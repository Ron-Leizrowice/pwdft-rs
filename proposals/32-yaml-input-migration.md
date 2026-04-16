# Proposal 32: Migrate to Pure YAML Input

## Problem

The codebase has two parallel input parsing paths that describe the same physical calculation:

1. **`src/input.rs`** — TOML-based `InputFile`, used by `main.rs` and all 4 example files (`examples/*.toml`)
2. **`src/settings.rs`** — YAML-based `Settings`, with 1 example file (`examples/si_scf_settings.yaml`) but **never used by `main.rs`**

The TOML path is the only one wired to the CLI. The YAML path is strictly better:

| Concern | TOML (`InputFile`) | YAML (`Settings`) |
|---|---|---|
| Sections | 3 (`system`, `kpoints`, `scf`) | 9 (`system`, `basis`, `kpoints`, `scf`, `electrons`, `xc`, `symmetry`, `pseudopotentials`, `output`) |
| Defaults | Partial — `ecut` lives in `[system]`, mixing in `[scf]` | Full — every optional section has `#[serde(default)]` with QE-convention defaults |
| Smearing types | None — only `smearing_sigma: f64` | Enum: `fermi_dirac`, `gaussian`, `methfessel_paxton`, `cold`, `fixed` |
| XC functional | Not configurable | Enum: `pz`, `pbe`, `pbe0`, `hse06` |
| Symmetry control | Not exposed | `enabled`, `time_reversal`, `tolerance` |
| Output control | Not exposed | `verbosity`, `write_density`, `write_bands` |
| Conversion helpers | `to_crystal()`, `to_high_sym_path()` | All of those plus `to_scf_params()`, `to_symmetry_info()`, `pseudopotential_path()`, `ecutwfc()`, `mp_grid()` |
| Test coverage | 1 test | 18 tests (roundtrip, defaults, partial overrides, error cases) |
| `n_bands` location | `system.n_bands` (physics leak) | `scf.n_bands` (correct section) |
| `ecut` location | `system.ecut` (mixed concerns) | `basis.ecutwfc` (own section) |
| Pseudopotentials | Nested under `[scf.pseudopotentials]` | Top-level `pseudopotentials:` section |

The TOML format also has ergonomic problems: TOML's `[[system.atoms]]` array-of-tables syntax is verbose and harder to read than YAML's list syntax for atom definitions.

Maintaining two parallel parsers doubles the surface area for bugs and means every new field must be added in two places. Only the weaker parser is actually used.

## Implementation

### Step 1: Wire `Settings` into `main.rs`

Replace `InputFile::from_file` with format detection based on file extension:

```rust
let config = match cli.input.extension().and_then(|e| e.to_str()) {
    Some("yaml" | "yml") => Settings::from_yaml_file(&cli.input)?,
    Some("toml") => Settings::from_toml_file(&cli.input)?,  // temporary bridge
    _ => return Err(PwdftError::Parse(
        "input file must have .yaml or .toml extension".into()
    )),
};
```

Refactor `main.rs` to use `Settings` as the single config type. The TOML branch would parse via `toml` then convert to `Settings` internally (temporary compatibility shim).

**Files:** `src/main.rs`, `src/settings.rs` (add `from_toml_file` bridge method)

### Step 2: Add TOML-to-Settings bridge

Add a `from_toml_str` / `from_toml_file` method on `Settings` that parses the old TOML format and maps it into the `Settings` struct. This provides backwards compatibility while we migrate example files:

```rust
impl Settings {
    pub fn from_toml_file(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        let old: InputFile = toml::from_str(&contents)
            .map_err(|e| PwdftError::Parse(e.to_string()))?;
        Ok(Self::from_legacy_input(old))
    }

    fn from_legacy_input(input: InputFile) -> Self {
        // Map InputFile fields → Settings sections
        // ...
    }
}
```

**Files:** `src/settings.rs`

### Step 3: Convert example files to YAML

Convert all 4 TOML examples to YAML equivalents:

| Old | New |
|---|---|
| `examples/si_scf.toml` | `examples/si_scf.yaml` |
| `examples/si_scf_converged.toml` | `examples/si_scf_converged.yaml` |
| `examples/si_scf_qe_match.toml` | `examples/si_scf_qe_match.yaml` |
| `examples/si_free_electron.toml` | `examples/si_free_electron.yaml` |

The existing `examples/si_scf_settings.yaml` serves as the template. Each new file should use the full sectioned layout with comments explaining each parameter.

**Files:** `examples/*.yaml` (4 new files), delete `examples/*.toml` (4 files)

### Step 4: Remove `input.rs` and TOML dependency

Once all examples and docs reference YAML:

1. Delete `src/input.rs`
2. Remove `pub mod input;` from `src/lib.rs`
3. Remove the `toml` crate from `Cargo.toml`
4. Remove the TOML bridge code added in Step 2
5. Update `Cli` struct doc comment from "Path to TOML input file" to "Path to YAML input file"
6. Update CLAUDE.md references (run commands, conventions section)

**Files:** `src/input.rs` (delete), `src/lib.rs`, `Cargo.toml`, `src/main.rs`, `CLAUDE.md`

### Step 5: Update documentation and CLI help

- Update `CLAUDE.md`: change all TOML references to YAML, update example commands
- Update `Cli` struct `#[arg]` help text
- Verify no stale `.toml` references remain anywhere in `src/` or `tests/`

**Files:** `CLAUDE.md`, `src/main.rs`

## Verification

1. **Unit tests pass:** `cargo test` — the 18 existing tests in `settings.rs` already cover YAML parsing, defaults, roundtrip, and error handling
2. **Example files parse:** Run `cargo run --release -- --input examples/si_scf.yaml` and verify SCF converges to the same energy as the TOML version
3. **Band structure works:** `cargo run --release -- --input examples/si_free_electron.yaml -o bands.tsv`
4. **No TOML remnants:** `grep -r '\.toml' src/ tests/` returns nothing (except Cargo.toml itself)
5. **Clippy clean:** `cargo clippy -q --all-targets`

## Estimated Effort

Straightforward — the hard work (designing the YAML schema, writing `Settings` with all sections, and 18 tests) is already done. This is mostly wiring and cleanup. A single focused session.
