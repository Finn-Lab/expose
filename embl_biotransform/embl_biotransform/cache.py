"""
cache.py
========

A persistent, process-safe cache on disk, so the pipeline stops depending on
the remote APIs being up at the moment you happen to ask.

Why this exists: ChEMBL's REST API has been returning 500s for long stretches,
MGnify analysis downloads take 5-15 seconds each, and a single comparison of
several analyses re-fetches all of it. Once something has been fetched
successfully it should never need fetching again.

Two layers use it:

- `install_http_cache()` hooks it into `fetchers._get_json` / `_get_text_file`,
  so *every* API response and results file is cached by URL. This is the layer
  that makes an outage survivable: whatever you fetched before an outage is
  still there during one.
- Callers cache composites directly (`get_or_set`), e.g. a built set of ChEMBL
  reference associations, which is many requests' worth of work.

Only successes are stored -- an error is never cached, so a failed call is
retried next time rather than remembered as failure.

Stored in SQLite (one file, no server, safe across threads and processes).
The default location follows the platform convention and can be overridden
with `$EXPOSE_CACHE_DIR`.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import sys
import threading
import time
from pathlib import Path
from typing import Any, Callable

#: Never expire by default: an MGnify analysis or a ChEMBL assay is immutable
#: once published. Searches pass a shorter ttl explicitly.
DEFAULT_TTL: float | None = None

SEARCH_TTL = 24 * 3600.0


def default_cache_dir() -> Path:
    """`$EXPOSE_CACHE_DIR`, else the platform's usual cache location."""
    override = os.environ.get("EXPOSE_CACHE_DIR")
    if override:
        return Path(override).expanduser()
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "expose"
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or (Path.home() / "AppData" / "Local")
        return Path(base) / "expose" / "cache"
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "expose"


class Cache:
    """Namespaced key/value store of JSON-serialisable values."""

    def __init__(self, path: str | Path | None = None, default_ttl: float | None = DEFAULT_TTL):
        if path is None:
            path = default_cache_dir() / "expose-cache.sqlite"
        self.path = Path(path)
        self.default_ttl = default_ttl
        self._local = threading.local()
        self._write_lock = threading.Lock()
        if str(self.path) != ":memory:":
            self.path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    # -- plumbing ------------------------------------------------------------ #

    @property
    def _db(self) -> sqlite3.Connection:
        """One connection per thread: sqlite3 connections aren't shareable."""
        conn = getattr(self._local, "conn", None)
        if conn is None:
            conn = sqlite3.connect(self.path, timeout=30)
            conn.execute("PRAGMA journal_mode=WAL")  # concurrent readers during a write
            conn.execute("PRAGMA busy_timeout=30000")
            self._local.conn = conn
        return conn

    def _init_db(self) -> None:
        with self._write_lock:
            self._db.execute("""
                CREATE TABLE IF NOT EXISTS entries (
                    namespace TEXT NOT NULL,
                    key       TEXT NOT NULL,
                    value     TEXT NOT NULL,
                    created   REAL NOT NULL,
                    expires   REAL,
                    PRIMARY KEY (namespace, key))""")
            self._db.commit()

    @staticmethod
    def make_key(*parts: Any) -> str:
        """A stable key for anything JSON-serialisable (URLs, params, id lists)."""
        blob = json.dumps(parts, sort_keys=True, default=str)
        return hashlib.sha256(blob.encode()).hexdigest()

    # -- api ------------------------------------------------------------------ #

    def get(self, namespace: str, key: str) -> Any | None:
        row = self._db.execute(
            "SELECT value, expires FROM entries WHERE namespace=? AND key=?",
            (namespace, key)).fetchone()
        if row is None:
            return None
        value, expires = row
        if expires is not None and expires < time.time():
            self.delete(namespace, key)
            return None
        try:
            return json.loads(value)
        except json.JSONDecodeError:  # corrupt row: drop it rather than fail
            self.delete(namespace, key)
            return None

    def set(self, namespace: str, key: str, value: Any, ttl: float | None = ...) -> None:
        ttl = self.default_ttl if ttl is ... else ttl
        expires = None if ttl is None else time.time() + ttl
        payload = json.dumps(value)
        with self._write_lock:
            self._db.execute(
                "INSERT OR REPLACE INTO entries (namespace, key, value, created, expires) "
                "VALUES (?,?,?,?,?)", (namespace, key, payload, time.time(), expires))
            self._db.commit()

    def get_or_set(self, namespace: str, key: str, fn: Callable[[], Any],
                   ttl: float | None = ...) -> Any:
        """Cached value, else call `fn` and store its result. An exception from
        `fn` propagates and nothing is stored."""
        hit = self.get(namespace, key)
        if hit is not None:
            return hit
        value = fn()
        if value is not None:
            self.set(namespace, key, value, ttl)
        return value

    def keys(self, namespace: str) -> list[str]:
        """Every live key in a namespace, newest first."""
        now = time.time()
        rows = self._db.execute(
            "SELECT key FROM entries WHERE namespace=? AND (expires IS NULL OR expires > ?) "
            "ORDER BY created DESC", (namespace, now)).fetchall()
        return [r[0] for r in rows]

    def items(self, namespace: str, limit: int | None = None) -> list[tuple[str, Any, float]]:
        """`(key, value, created)` for a namespace, newest first.

        Only usable where keys are meaningful (the ChEMBL library keys entries
        by assay accession); namespaces keyed by `make_key` hashes cannot be
        listed back into anything readable.
        """
        now = time.time()
        sql = ("SELECT key, value, created FROM entries WHERE namespace=? "
               "AND (expires IS NULL OR expires > ?) ORDER BY created DESC")
        params: tuple = (namespace, now)
        if limit:
            sql += " LIMIT ?"
            params = (namespace, now, limit)
        out = []
        for key, value, created in self._db.execute(sql, params).fetchall():
            try:
                out.append((key, json.loads(value), created))
            except json.JSONDecodeError:
                continue
        return out

    def delete(self, namespace: str, key: str) -> None:
        with self._write_lock:
            self._db.execute("DELETE FROM entries WHERE namespace=? AND key=?", (namespace, key))
            self._db.commit()

    def clear(self, namespace: str | None = None) -> int:
        with self._write_lock:
            cur = (self._db.execute("DELETE FROM entries WHERE namespace=?", (namespace,))
                   if namespace else self._db.execute("DELETE FROM entries"))
            self._db.commit()
            return cur.rowcount

    def purge_expired(self) -> int:
        with self._write_lock:
            cur = self._db.execute("DELETE FROM entries WHERE expires IS NOT NULL AND expires < ?",
                                   (time.time(),))
            self._db.commit()
            return cur.rowcount

    def stats(self) -> dict:
        rows = self._db.execute(
            "SELECT namespace, COUNT(*), SUM(LENGTH(value)) FROM entries GROUP BY namespace").fetchall()
        size = 0
        if str(self.path) != ":memory:" and self.path.exists():
            size = self.path.stat().st_size
        return {
            "path": str(self.path),
            "file_bytes": size,
            "entries": sum(r[1] for r in rows),
            "namespaces": {r[0]: {"entries": r[1], "bytes": r[2] or 0} for r in rows},
        }

    def close(self) -> None:
        conn = getattr(self._local, "conn", None)
        if conn is not None:
            conn.close()
            self._local.conn = None


# --------------------------------------------------------------------------- #
# HTTP layer
# --------------------------------------------------------------------------- #

def install_http_cache(cache: Cache) -> None:
    """Route `fetchers._get_json` / `_get_text_file` through `cache`.

    Wrapping the transport rather than each client means every source (MGnify,
    ChEMBL, InterPro, the FTP result files) is covered, including calls made
    deep inside a build.
    """
    from . import fetchers

    if getattr(fetchers, "_cache_installed", False):
        fetchers.HTTP_CACHE = cache
        return

    raw_json, raw_text = fetchers._get_json, fetchers._get_text_file

    def cached_json(url: str, params: dict | None = None, **kwargs):
        active = fetchers.HTTP_CACHE
        if active is None:
            return raw_json(url, params=params, **kwargs)
        ttl = SEARCH_TTL if _is_search(url, params) else ...
        return active.get_or_set("http.json", Cache.make_key(url, params),
                                 lambda: raw_json(url, params=params, **kwargs), ttl)

    def cached_text(url: str, **kwargs):
        active = fetchers.HTTP_CACHE
        if active is None:
            return raw_text(url, **kwargs)
        return active.get_or_set("http.text", Cache.make_key(url),
                                 lambda: raw_text(url, **kwargs))

    fetchers._get_json, fetchers._get_text_file = cached_json, cached_text
    fetchers.HTTP_CACHE = cache
    fetchers._cache_installed = True


def _is_search(url: str, params: dict | None) -> bool:
    """Searches and listings can gain new hits, so they expire; records keyed by
    accession never do."""
    if params and any(k in params for k in ("search", "title", "assay_organism__icontains",
                                            "assay_description__icontains")):
        return True
    return url.rstrip("/").endswith(("/studies", "/analyses", "/samples"))
