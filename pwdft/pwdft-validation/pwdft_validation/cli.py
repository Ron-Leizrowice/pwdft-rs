"""Dispatcher for the pwdft-validate CLI.

Usage:
    pwdft-validate <script-name> [script args...]
    pwdft-validate --list
    pwdft-validate --help

``<script-name>`` matches a module under ``pwdft_validation.scripts`` (e.g.
``dij_reference``, ``vgch_sad_heavy``). The dispatcher sets ``sys.argv`` and
runs the module's ``main`` callable.

Existing scripts are preserved as-is and can also be run directly:

    uv run python -m pwdft_validation.scripts.dij_reference --help
"""

from __future__ import annotations

import argparse
import importlib
import pkgutil
import sys
from collections.abc import Sequence

from pwdft_validation import scripts as scripts_pkg


def _available_scripts() -> list[str]:
    return sorted(mod.name for mod in pkgutil.iter_modules(scripts_pkg.__path__) if not mod.name.startswith("_"))


def _run_script(name: str, args: Sequence[str]) -> int:
    full = f"pwdft_validation.scripts.{name}"
    try:
        mod = importlib.import_module(full)
    except ImportError as e:
        print(f"pwdft-validate: no such script '{name}' ({e})", file=sys.stderr)
        return 2

    main_fn = getattr(mod, "main", None)
    if main_fn is None:
        print(
            f"pwdft-validate: script '{name}' has no main() entry point",
            file=sys.stderr,
        )
        return 2

    sys.argv = [f"pwdft-validate {name}", *args]
    rc = main_fn()
    if rc is None:
        return 0
    if isinstance(rc, int):
        return rc
    print(f"pwdft-validate: script '{name}' returned non-int {rc!r}", file=sys.stderr)
    return 1


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pwdft-validate",
        description="Run validation / reference-data scripts for pwdft-rs.",
        allow_abbrev=False,
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available scripts and exit.",
    )
    parser.add_argument(
        "script",
        nargs="?",
        help="Name of the script module under pwdft_validation.scripts.",
    )
    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to the script.",
    )
    args = parser.parse_args(argv)

    if args.list:
        for name in _available_scripts():
            print(name)
        return 0

    if args.script is None:
        parser.print_help()
        return 0

    return _run_script(args.script, args.script_args)


if __name__ == "__main__":
    raise SystemExit(main())
