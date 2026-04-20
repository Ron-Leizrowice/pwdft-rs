# pwdft-validation

Quantum ESPRESSO cross-checks and reference-data generators for
[pwdft-rs](../README.md).

The workspace's Rust test suite compares pwdft-rs results against
QE reference outputs stored here. This package owns:

- The **reference data** those Rust tests consume
  (`reference/csv/*.csv`, `reference/qe/*.{in,out,bin,toml}`).
- The **generators** that produced it — Python scripts under
  `src/pwdft_validation/scripts/`. Each script either parses a QE
  output file or runs pwdft-rs and writes a CSV pin.
- A shared `pwdft_validation` package with filesystem constants
  (`REPO_ROOT`, `QE_REF_DIR`, `CSV_REF_DIR`, …) so scripts don't
  hand-roll `Path(__file__).parents[N]` chains.

## Layout

```
validation/
├── pyproject.toml              ← uv / hatchling project
├── README.md
├── src/pwdft_validation/
│   ├── __init__.py             ← re-exports paths module
│   ├── _paths.py               ← REPO_ROOT, QE_REF_DIR, …
│   ├── cli.py                  ← pwdft-validate dispatcher
│   └── scripts/
│       ├── dij_reference.py
│       ├── vgch_sad_heavy.py
│       └── …                   ← 15+ generators
├── reference/
│   ├── csv/                    ← Rust tests read these
│   └── qe/                     ← QE pw.x inputs + outputs + binaries
└── tests/                      ← Python unit tests
```

## Running

```bash
# Install in editable mode (dev extras: ruff, mypy, pytest).
cd validation/
uv sync --extra dev

# Run a named script (forwards remaining args):
uv run pwdft-validate vgch_sad_heavy --help

# Or directly:
uv run python -m pwdft_validation.scripts.vgch_sad_heavy

# List available scripts:
uv run pwdft-validate --list
```

## Developing

```bash
uv run ruff check .
uv run ruff format .
uv run mypy src tests
uv run pytest
```

CI (GitHub Actions) runs all four on pull requests — see
`.github/workflows/validation.yml`.

## Relationship to the Rust suite

Rust tests locate reference data via `env!("CARGO_WORKSPACE_DIR")`:

```rust
let csv = PathBuf::from(env!("CARGO_WORKSPACE_DIR"))
    .join("validation/reference/csv/foo.csv");
```

When regenerating a CSV pin, run the corresponding generator script
from this package, then re-run the Tier-1 Rust test that consumes it.
