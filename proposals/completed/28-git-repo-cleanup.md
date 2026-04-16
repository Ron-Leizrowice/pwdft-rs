# Proposal 28: Git Repository Cleanup

## Problem

The git repository has accumulated several issues that bloat its size and leave broken references:

### 1. Partial download files in git history (~19 MB)

Two `pseudopotentials/Unconfirmed*` files (10 MB + 9 MB) — partial browser downloads — were accidentally committed in 568f224 and removed in ac98a65. They are no longer tracked but remain in the git pack, permanently inflating clone size.

### 2. USPP and PAW pseudopotentials can't be used (48 MB tracked)

The code only supports norm-conserving (NC) pseudopotentials via the UPF parser. USPP and PAW support does not exist and is not planned near-term. Yet 53 USPP/PAW files (48 MB) are tracked:

| Directory | Files | Size | Used by code? |
|-----------|-------|------|--------------|
| `pseudopotentials/nc/lda/` | 74 | 17 MB | **Yes** (Si, Fe, C used in tests) |
| `pseudopotentials/nc/pbe/` | 72 | 18 MB | **No** (PBE XC not implemented) |
| `pseudopotentials/uspp/pbe/` | 41 | 28 MB | **No** (USPP not supported) |
| `pseudopotentials/paw/pbe/` | 12 | 20 MB | **No** (PAW not supported) |

Only 3 files are actually used: `nc/lda/Si.upf`, `nc/lda/Fe.upf`, `nc/lda/C.upf`.

### 3. Broken `pseudopotentials/Si.UPF` references

8 files reference `pseudopotentials/Si.UPF` which no longer exists (removed when PPs were reorganized into `nc/lda/`). This breaks benchmarks and examples:

| File | Line |
|------|------|
| `benches/scf_benchmarks.rs` | 44 |
| `examples/si_scf.toml` | 34 |
| `examples/si_scf_qe_match.toml` | 34 |
| `examples/si_scf_converged.toml` | 32 |
| `examples/si_scf_settings.yaml` | 46 |
| `src/settings.rs` | 496, 555, 629 (test YAML fixtures) |

### 4. PBE pseudopotentials without PBE XC

72 PBE norm-conserving PPs (18 MB) are tracked but the code only implements LDA exchange-correlation. Using an LDA XC functional with a PBE pseudopotential is technically valid but suboptimal — and no test or example uses these files.

## Implementation

### Step 1: Fix broken `Si.UPF` references

Replace all `pseudopotentials/Si.UPF` paths with `pseudopotentials/nc/lda/Si.upf`:

```
benches/scf_benchmarks.rs:  "pseudopotentials/Si.UPF" → "pseudopotentials/nc/lda/Si.upf"
examples/si_scf.toml:       Si = "../pseudopotentials/Si.UPF" → Si = "pseudopotentials/nc/lda/Si.upf"
examples/si_scf_qe_match.toml: same
examples/si_scf_converged.toml: same
examples/si_scf_settings.yaml: same
src/settings.rs:            test YAML fixtures (3 occurrences)
```

Note: example TOML files use relative paths from the example directory (`../pseudopotentials/`). After the fix, they should use paths relative to the project root, matching how `cargo run` resolves them.

### Step 2: Remove unused pseudopotentials

Remove from tracking (not from git history — that requires rewriting):

```bash
git rm -r pseudopotentials/uspp/
git rm -r pseudopotentials/paw/
git rm -r pseudopotentials/nc/pbe/
```

This removes 125 files (~66 MB) that the code cannot use. Keep only `pseudopotentials/nc/lda/` (74 files, 17 MB) which contains the PPs actually used in tests and referenced by examples.

Update `pseudopotentials/README.md` to remove SSSP/USPP/PAW sections and note that PBE PPs can be re-added when GGA is implemented.

### Step 3: Trim unused NC/LDA pseudopotentials (optional)

Of the 74 `nc/lda/` files, only 3 are used in tests (Si, Fe, C). The rest (71 files, ~16 MB) are there for user convenience. Options:

**Option A (conservative):** Keep all 74. Users may want other elements without downloading separately.

**Option B (minimal):** Keep only the 3-4 used by tests plus a handful of common elements (Al, O, N, Cu, Ga — used in typical DFT tutorials). Document download instructions for the full PseudoDojo set.

Recommend **Option A** for now — NC/LDA PPs are the ones the code actually supports.

### Step 4: Purge partial downloads from history (optional, destructive)

The `Unconfirmed*` files add ~19 MB to every clone. To remove them:

```bash
git filter-repo --invert-paths --path-glob 'pseudopotentials/Unconfirmed*'
```

This rewrites history. Only do this if:
- No other clones/forks exist that would be broken
- You're comfortable with force-pushing to origin

If you prefer not to rewrite history, the files will persist in the pack but have no effect on the working tree. `.gitignore` already prevents new ones via `*.tar.gz` / `*.tgz` — add `Unconfirmed*` as well.

### Step 5: Update `.gitignore`

Add patterns to prevent future accidents:

```gitignore
# Partial browser downloads
Unconfirmed*
*.crdownload
*.part

# Pseudopotential formats we don't support yet
pseudopotentials/uspp/
pseudopotentials/paw/
pseudopotentials/nc/pbe/
```

## Impact

| Change | Size saved | Risk |
|--------|-----------|------|
| Fix Si.UPF references | 0 | Zero — fixes broken paths |
| Remove USPP/PAW/PBE PPs | ~66 MB from working tree | Zero — code can't use them |
| Purge Unconfirmed from history | ~19 MB from pack | Medium — requires history rewrite |
| Update .gitignore | 0 | Zero |

## Acceptance Criteria

1. `cargo bench` compiles (currently broken due to missing `Si.UPF`).
2. Example TOML files point to existing PP paths.
3. `pseudopotentials/uspp/` and `pseudopotentials/paw/` no longer tracked.
4. `pseudopotentials/nc/pbe/` no longer tracked.
5. `.gitignore` prevents `Unconfirmed*` and partial downloads.
6. `cargo test` passes — no PP path regressions.
7. Repository working tree drops from ~82 MB to ~17 MB in pseudopotentials.
