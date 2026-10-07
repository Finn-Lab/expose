"""
expose_app
==========

The EXPOSE application shell: the CLI, and the native desktop window.

The analysis library (`embl_biotransform`) and the HTTP API (`webapp`) know
nothing about this layer, so the same server runs behind a desktop window, in
a browser, or headless on a server.
"""

__all__ = ["cli", "desktop"]
