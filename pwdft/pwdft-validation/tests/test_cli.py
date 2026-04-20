"""Smoke tests for the pwdft-validate cyclopts CLI."""

from __future__ import annotations

import pytest
from pwdft_validation import CSV_REF_DIR, PROJECT_ROOT, PSEUDO_DIR, QE_REF_DIR
from pwdft_validation.cli import app
from pwdft_validation.density import standalone_app as density_standalone_app
from pwdft_validation.diagnostics.residual_scan import standalone_app as residual_scan_standalone_app


def test_paths_point_at_workspace_dirs() -> None:
    assert PROJECT_ROOT.is_dir()
    assert (PROJECT_ROOT / "Cargo.toml").is_file()
    assert QE_REF_DIR.is_dir()
    assert CSV_REF_DIR.is_dir()
    assert PSEUDO_DIR.is_dir()


def test_app_help_exits_zero(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        app(["--help"])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "pwdft-validate" in captured.out


def test_app_registers_reference_subcommands() -> None:
    """Every top-level sub-app and command must be registered on the root app."""
    registered = set(app._commands.keys()) if hasattr(app, "_commands") else set()
    # Fall back to iterating the meta / help text if the private API changes.
    if not registered:
        with pytest.raises(SystemExit):
            app(["--help"])
        return
    expected = {"reference", "energy", "diag", "fermi", "density", "pbe"}
    assert expected.issubset(registered), f"missing subcommands: {expected - registered}"


def test_reference_subapp_help(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        app(["reference", "--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    for cmd in ("vloc", "beta", "dij", "hamiltonian", "sad", "nlcc"):
        assert cmd in out, f"sub-command {cmd!r} missing from `reference --help`"


def test_density_standalone_help(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        density_standalone_app(["--help"])
    assert exc.value.code == 0
    assert "pwdft-density-parse" in capsys.readouterr().out


def test_residual_scan_standalone_help(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        residual_scan_standalone_app(["--help"])
    assert exc.value.code == 0
    assert "pwdft-residual-scan" in capsys.readouterr().out


def test_unknown_command_errors() -> None:
    with pytest.raises(SystemExit) as exc:
        app(["nonexistent_subcommand_xyz"])
    assert exc.value.code != 0
