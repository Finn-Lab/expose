"""
jobs.py
=======

A small in-process job queue, so slow or flaky work stops blocking the UI.

The web viewer used to call ChEMBL and MGnify synchronously inside a request.
A ChEMBL search can take two filter queries with three retries each, and an
MGnify analysis download 5-15 seconds; during a ChEMBL outage a request could
hang for minutes and then fail, with nothing to show for it. Now the request
submits a job and returns immediately, the front end polls for progress, and
partial results survive because each step writes through the cache.

Deliberately modest: threads, in one process, with no persistence of the queue
itself. What *is* persistent is everything a job fetches, via `cache.py` -- so
a job that dies halfway leaves its completed fetches behind and a re-run picks
up from there. Jobs are for progress and responsiveness; the cache is for
durability.
"""

from __future__ import annotations

import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Callable

PENDING, RUNNING, DONE, FAILED, CANCELLED = "pending", "running", "done", "failed", "cancelled"
TERMINAL = (DONE, FAILED, CANCELLED)


class JobCancelled(Exception):
    """Raised inside a job when the caller has asked it to stop."""


@dataclass
class Job:
    id: str
    name: str
    state: str = PENDING
    progress: float = 0.0          # 0..1, -1 when the total isn't known
    message: str = ""
    result: Any = None
    error: str | None = None
    meta: dict = field(default_factory=dict)
    submitted: float = field(default_factory=time.time)
    started: float | None = None
    finished: float | None = None
    _cancel: threading.Event = field(default_factory=threading.Event, repr=False)

    def to_dict(self, include_result: bool = True) -> dict:
        out = {
            "id": self.id, "name": self.name, "state": self.state,
            "progress": self.progress, "message": self.message,
            "error": self.error, "meta": self.meta,
            "submitted": self.submitted, "started": self.started, "finished": self.finished,
            "elapsed": round((self.finished or time.time()) - (self.started or self.submitted), 2),
        }
        if include_result and self.state == DONE:
            out["result"] = self.result
        return out


class JobContext:
    """Handed to the job function: report progress, and check for cancellation."""

    def __init__(self, job: Job):
        self._job = job

    def progress(self, fraction: float | None = None, message: str | None = None) -> None:
        self.check_cancelled()
        if fraction is not None:
            self._job.progress = max(0.0, min(1.0, fraction))
        if message is not None:
            self._job.message = message

    def step(self, done: int, total: int, message: str | None = None) -> None:
        self.progress(done / total if total else -1.0, message)

    def check_cancelled(self) -> None:
        if self._job._cancel.is_set():
            raise JobCancelled()

    @property
    def cancelled(self) -> bool:
        return self._job._cancel.is_set()


class JobQueue:
    """Submit callables, poll their state. Finished jobs are kept for a while
    so a page reload can still collect a result it missed."""

    def __init__(self, max_workers: int = 4, keep_finished: int = 100,
                 retain_seconds: float = 3600.0):
        self._pool = ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="expose-job")
        self._jobs: dict[str, Job] = {}
        self._order: list[str] = []
        self._lock = threading.Lock()
        self.keep_finished = keep_finished
        self.retain_seconds = retain_seconds

    def submit(self, name: str, fn: Callable[..., Any], *args, meta: dict | None = None,
               **kwargs) -> Job:
        """`fn` is called as `fn(ctx, *args, **kwargs)` if it accepts a context
        as its first argument -- pass `wants_context=False` in `meta` if not."""
        job = Job(id=uuid.uuid4().hex[:12], name=name, meta=dict(meta or {}))
        with self._lock:
            self._jobs[job.id] = job
            self._order.append(job.id)
            self._prune_locked()
        self._pool.submit(self._run, job, fn, args, kwargs)
        return job

    def _run(self, job: Job, fn, args, kwargs) -> None:
        if job._cancel.is_set():
            job.finished, job.state = time.time(), CANCELLED
            return
        job.state, job.started = RUNNING, time.time()
        ctx = JobContext(job)
        state, error = DONE, None
        try:
            job.result = fn(ctx, *args, **kwargs)
            job.progress = 1.0
        except JobCancelled:
            state = CANCELLED
            job.message = "cancelled"
        except Exception as exc:  # noqa: BLE001 - a job must never kill the worker
            state, error = FAILED, f"{type(exc).__name__}: {exc}"
            job.meta.setdefault("traceback", traceback.format_exc(limit=5))
        finally:
            # `finished` must land before `state`: another thread pruning the
            # queue reads state first, and a job that looks terminal with no
            # timestamp would be mistaken for an ancient one and dropped.
            job.error = error
            job.finished = time.time()
            job.state = state

    def get(self, job_id: str) -> Job | None:
        return self._jobs.get(job_id)

    def cancel(self, job_id: str) -> bool:
        """Ask a job to stop. It stops at its next progress report, so a job
        blocked in a long HTTP retry ends when that call returns."""
        job = self._jobs.get(job_id)
        if job is None or job.state in TERMINAL:
            return False
        job._cancel.set()
        if job.state == PENDING:
            job.state, job.finished = CANCELLED, time.time()
        return True

    def list(self, limit: int = 50) -> list[dict]:
        with self._lock:
            ids = self._order[-limit:][::-1]
        return [self._jobs[i].to_dict(include_result=False) for i in ids if i in self._jobs]

    def _prune_locked(self) -> None:
        """Forget old finished jobs; never touch anything still running."""
        now = time.time()
        cutoff = now - self.retain_seconds
        finished = [i for i in self._order
                    if i in self._jobs and self._jobs[i].state in TERMINAL]
        # a missing timestamp means "just now", never "long ago"
        stale = [i for i in finished if (self._jobs[i].finished or now) < cutoff]
        excess = finished[:-self.keep_finished] if len(finished) > self.keep_finished else []
        for job_id in set(stale) | set(excess):
            self._jobs.pop(job_id, None)
            if job_id in self._order:
                self._order.remove(job_id)

    def shutdown(self, wait: bool = False) -> None:
        for job in self._jobs.values():
            if job.state not in TERMINAL:
                job._cancel.set()
        self._pool.shutdown(wait=wait)
