"""`python -m webapp` is kept as an alias for the `expose` command, so older
instructions keep working; new entry points live in `expose_app.cli`."""

from expose_app.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
