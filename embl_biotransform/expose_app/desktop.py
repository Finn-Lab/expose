"""
desktop.py
==========

The native application window.

The Flask server runs on a daemon thread and pywebview owns the main thread,
which is a hard requirement on macOS: the Cocoa event loop must be on thread
one. Closing the window ends the process, and the daemon thread goes with it.

pywebview is an optional dependency. Importing this module without it raises
ImportError, which `cli.py` catches and falls back to a browser tab -- so the
app still runs everywhere, just without its own window.
"""

from __future__ import annotations

import threading
import time
import urllib.error
import urllib.request

import webview  # noqa: F401  -- ImportError here is the signal to fall back

WINDOW_TITLE = "EXPOSE - biotransformation evidence viewer"
MIN_SIZE = (1100, 720)


def _serve(app, host: str, port: int, debug: bool) -> None:
    app.run(host=host, port=port, debug=debug, use_reloader=False, threaded=True)


def _wait_until_up(url: str, timeout: float = 15.0) -> bool:
    """Don't show the window until the server answers, or it opens on an error
    page and the user has to reload."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(f"{url}/api/health", timeout=1):
                return True
        except urllib.error.HTTPError:
            return True          # responding, even if unhappy
        except Exception:        # noqa: BLE001 - not up yet
            time.sleep(0.15)
    return False


def run_window(app, host: str, port: int, url: str, debug: bool = False) -> int:
    server = threading.Thread(target=_serve, args=(app, host, port, debug),
                              name="expose-server", daemon=True)
    server.start()

    if not _wait_until_up(url):
        print("  server did not start in time; not opening the window")
        return 1

    webview.create_window(WINDOW_TITLE, url,
                          width=max(MIN_SIZE[0], 1280), height=max(MIN_SIZE[1], 860),
                          min_size=MIN_SIZE, text_select=True)
    webview.start()              # blocks on the main thread until the window closes
    return 0
