"""Smoke tests for the pwdft-validate CLI dispatcher."""

from __future__ import annotations

import pytest
from pwdft_validation import CSV_REF_DIR, PSEUDO_DIR, QE_REF_DIR, REPO_ROOT
from pwdft_validation.cli import _available_scripts, main


def test_paths_point_at_workspace_dirs() -> None:
    assert REPO_ROOT.is_dir()
    assert (REPO_ROOT / "Cargo.toml").is_file()
    assert QE_REF_DIR.is_dir()
    assert CSV_REF_DIR.is_dir()
    assert PSEUDO_DIR.is_dir()


def test_list_nonempty() -> None:
    scripts = _available_scripts()
    assert scripts, "expected at least one validation script to be discoverable"
    # Private/dunder modules are excluded.
    assert not any(s.startswith("_") for s in scripts)


def test_help_exits_zero(capsys: pytest.CaptureFixture[str]) -> None:
    rc = main([])
    assert rc == 0
    captured = capsys.readouterr()
    assert "pwdft-validate" in captured.out


def test_list_flag(capsys: pytest.CaptureFixture[str]) -> None:
    rc = main(["--list"])
    assert rc == 0
    listed = capsys.readouterr().out.strip().splitlines()
    assert listed, "expected --list to print at least one script name"


def test_unknown_script_returns_nonzero(capsys: pytest.CaptureFixture[str]) -> None:
    rc = main(["nonexistent_script_name"])
    assert rc != 0
    assert "no such script" in capsys.readouterr().err
