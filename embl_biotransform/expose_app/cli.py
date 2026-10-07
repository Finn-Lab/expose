"""
cli.py
======

The `expose` command.

    expose                    # desktop window (falls back to a browser tab)
    expose --browser          # force a browser tab
    expose --headless         # serve only; for a remote or shared machine
    expose cache stats        # what has been cached, and where
    expose cache clear [ns]   # drop cached responses

The server binds to 127.0.0.1 by default: this is a local analysis tool with
no authentication, and the cache holds whatever you have fetched.
"""

from __future__ import annotations

import argparse
import sys
import threading
import webbrowser

from embl_biotransform.cache import Cache, default_cache_dir, install_http_cache


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="expose", description="EXPOSE - biotransformation evidence viewer")
    parser.add_argument("--host", default="127.0.0.1",
                        help="interface to bind (default 127.0.0.1; use 0.0.0.0 to share)")
    parser.add_argument("--port", type=int, default=8050)
    parser.add_argument("--browser", action="store_true", help="open a browser tab instead of the app window")
    parser.add_argument("--headless", action="store_true", help="serve only: no window, no browser")
    parser.add_argument("--no-cache", action="store_true", help="do not read or write the on-disk cache")
    parser.add_argument("--cache-dir", help="override the cache location")
    parser.add_argument("--workers", type=int, default=4, help="background job workers (default 4)")
    parser.add_argument("--debug", action="store_true")

    sub = parser.add_subparsers(dest="command")
    cache_cmd = sub.add_parser("cache", help="inspect or clear the on-disk cache")
    cache_cmd.add_argument("action", choices=["stats", "clear", "purge", "path"])
    cache_cmd.add_argument("namespace", nargs="?", help="limit 'clear' to one namespace")
    return parser


def _open_cache(args) -> Cache | None:
    if getattr(args, "no_cache", False):
        return None
    path = None
    if getattr(args, "cache_dir", None):
        from pathlib import Path
        path = Path(args.cache_dir).expanduser() / "expose-cache.sqlite"
    return Cache(path)


def _cache_command(args) -> int:
    cache = _open_cache(argparse.Namespace(no_cache=False, cache_dir=args.cache_dir))
    if args.action == "path":
        print(cache.path)
    elif args.action == "stats":
        stats = cache.stats()
        print(f"cache: {stats['path']}")
        print(f"  {stats['entries']} entries, {stats['file_bytes'] / 1e6:.1f} MB on disk")
        for name, info in sorted(stats["namespaces"].items()):
            print(f"    {name:16} {info['entries']:>6} entries  {info['bytes'] / 1e6:>7.1f} MB")
        if not stats["namespaces"]:
            print("    (empty)")
    elif args.action == "clear":
        removed = cache.clear(args.namespace)
        print(f"removed {removed} entries" + (f" from {args.namespace}" if args.namespace else ""))
    elif args.action == "purge":
        print(f"removed {cache.purge_expired()} expired entries")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command == "cache":
        return _cache_command(args)

    cache = _open_cache(args)
    if cache is not None:
        install_http_cache(cache)

    from webapp.app import create_app
    app = create_app(cache=cache, max_workers=args.workers)

    url = f"http://{args.host}:{args.port}"
    print(f"EXPOSE running at {url}")
    print(f"  cache: {cache.path if cache else 'disabled'}")
    print("  Ctrl+C to stop")

    if args.headless:
        app.run(host=args.host, port=args.port, debug=args.debug, use_reloader=False)
        return 0

    if not args.browser:
        try:
            from .desktop import run_window
        except ImportError:
            print("  (pywebview not installed - opening a browser tab instead;"
                  " `pip install 'expose-biotransform[desktop]'` for the app window)")
        else:
            return run_window(app, host=args.host, port=args.port, url=url, debug=args.debug)

    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    app.run(host=args.host, port=args.port, debug=args.debug, use_reloader=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
