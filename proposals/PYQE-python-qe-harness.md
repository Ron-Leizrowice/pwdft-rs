---
id: PYQE
status: active
priority: high
complexity: large
risk: medium
depends_on: []
blocks: [VQEF, VGCH-2, VGCH-MECH]
---

# PYQE: Route all QE invocation through the Python validation module + CI regression gate

## Problem

QE is invoked from three disconnected places in the repo, each with its own orchestration contract. The resulting workflow is ad-hoc and brittle, and there is no automated gate against QE-reference regressions.

1. **Humans (and agents) shell out by hand.** `.claude/skills/qe-runner/SKILL.md` instructs every agent to hand-assemble an `mpirun pw.x` pipeline, wrap it in `.claude/bin/machine-lock run ...`, plumb `gtimeout`, detect cores per-OS, etc. Every QE run is a fresh copy-paste of a 3-line bash incantation. Nothing validates the output after the run — the agent eyeballs the `.out` file.
2. **Python parses but never runs.** `pwdft_validation.diagnostics.pbe_refs.extract` (and the VGCH-2 energy/fermi modules) only *read* pre-generated `data/qe/*.out` files. If the outputs are stale, the Python never notices — it re-parses yesterday's reference with no version check.
3. **Rust tests read a static TOML.** `pwdft/pwdft-core/tests/qe_validation.rs` consumes `data/qe/reference_data.toml`. That TOML was last hand-regenerated on 2026-04-17 + 2026-04-19 (see file header). There is no automation to re-regenerate on a QE version bump, pseudopotential change, or input-deck edit; regressions of QE-side data go undetected until an agent happens to re-run by hand.
4. **Lock policy is enforced by convention.** QELK (completed) made the rule "every QE run holds the machine lock" a CLAUDE.md policy. There is no code path that enforces it — an agent that forgets the wrapper silently corrupts benchmark wall-time measurements.
5. **CI ignores QE entirely.** `.github/workflows/ci.yml` runs clippy + `cargo test -q -p pwdft-core` (Tier 1 only) + ruff + ty + pytest. No job runs QE, no job re-validates the TOML pins, and no job gates a PR on VQEF residual regression. The 16-cell VQEF matrix (`4 GREEN / 12 YELLOW`) is scored by hand during grooming passes.

**Concrete evidence.**

| Source | QE entry point | Lock enforcement | Post-run validation |
|---|---|---|---|
| `.claude/skills/qe-runner/SKILL.md` | hand-typed `mpirun pw.x` | documented; human follows | none |
| `pwdft_validation.diagnostics.pbe_refs` | (reads cached `.out`) | n/a | regex scrape |
| `pwdft/pwdft-core/tests/qe_validation.rs` | (reads `reference_data.toml`) | n/a | assertions vs pwdft |
| CI (`.github/workflows/ci.yml`) | — | — | — |

`rg -l` confirms six files reference `mpirun`/`pw.x -in` directly; none share infrastructure.

**Why it matters now.** VQEF is the outstanding high-priority item (12 YELLOW cells, gating VGCH-MECH and Class A/B/C diagnostics). Every time we tighten a tolerance, re-run a reference, or add a new cell, we add another hand-coded shell snippet with its own bugs. RWHK C1 already flagged one case where the GPU-vs-CPU consistency test was structurally wrong because nothing automated cross-checked the invocation. A Python-orchestrated runner closes the same class of hole on the QE side.

## Research

### Current Python module layout (post-refactor — see commits `9c2157d`, `a19c0cd`)

```text
pwdft/pwdft-validation/pwdft_validation/
├── cli.py              cyclopts app — `pwdft-validate <sub> ...`
├── paths.py            PROJECT_ROOT / PSEUDO_DIR / QE_REF_DIR / CSV_REF_DIR
├── units.py            BOHR_TO_ANG, RY_TO_EV, ...
├── upf.py              shared UPF v2 parser
├── integrate.py        simpson_qe quadrature
├── energy.py           extract_components / extract_trace / join_traces
├── fermi.py            QE Fermi-finder reference
├── density.py          charge-density.dat → VGCH2BIN parser
├── reference/          vloc / beta / dij / hamiltonian / sad / nlcc generators
└── diagnostics/        ef_shift / residual_scan / pbe_refs
```

The module already owns the **parse/extract** side of the QE boundary. Adding a **run** sibling package is the natural next step.

### Lock-orchestration options

| Approach | Pros | Cons |
|---|---|---|
| Python calls `subprocess.run(["machine-lock", "run", ...])` | zero-new-code; reuses MLFX's atomic mkdir | leaks shell concerns into Python; pytest gets a shell dep |
| Python re-implements `mkdir`-atomic lock in `pwdft_validation.lock` | single language; easy to unit-test | duplicates MLFX's hardening (PID liveness, 3h cap, worktree owner check) |
| Python calls a thin Rust CLI in `pwdft-validation` crate | type-safe; shared with Rust tests | new crate, build-time dep for Python |

Recommendation: **Option 1 plus a thin `pwdft_validation.lock` context-manager wrapper** that invokes the existing shell script. Keeps a single source of truth for the atomic primitive (MLFX shell), gives Python a clean `with lock(...)` call site, preserves the cross-worktree hook enforcement. We can promote to Option 2 later if pytest-side unit tests demand a mock.

### QE in CI — build vs binary-distribution

GitHub Actions Ubuntu runners do not have QE preinstalled. Options (measured on `ubuntu-latest`, 4 vCPU):

| Strategy | First-run time | Cached time | Fragility |
|---|---|---|---|
| Build QE 7.5 from source (`make pw`, `-j4`) | ~18 min | ~2 min (artifact cache) | medium — libxc/fftw pinning |
| Install via `conda-forge::quantum-espresso` | ~2 min | ~40 s (conda cache) | low — stable recipe |
| Pull a prebuilt Docker image | ~1 min | ~15 s | low — image pinning |

Recommendation: **`conda-forge::quantum-espresso`** — fastest cached path, low maintenance, broadly trusted recipe. Build-from-source is available as a fallback job if conda drifts from upstream 7.5.

### Scope — which invocation patterns must PYQE cover

From the reference inputs + SKILL.md:

- **SCF (`pw.x`)** — 16 cells (Si/C/Al/Fe/Cu/GaAs/NaCl/MgO × LDA/PBE). Used by VQEF, BSUM, VGCH-2.
- **Post-processing** — `pp.x` (charge densities — `charge-density.dat` → VGCH2BIN transplant), `bands.x` (band structure), `projwfc.x` (DOS / PDOS).
- **Phonons** — `ph.x`, `q2r.x`, `matdyn.x`. Not yet used; bundle anyway so the skill covers it.

PYQE's v1 delivers SCF + charge-density. Bands/PDOS/phonons land under the same harness on demand.

## Implementation

### Phase A — Python QE runner package (`pwdft_validation.qe`)

Add a new sub-package `pwdft_validation/qe/` with:

- `pwdft_validation/qe/__init__.py` — re-exports the public API.
- `pwdft_validation/qe/binary.py` — locates `pw.x`/`pp.x`/etc. Tries env overrides first (`QE_PREFIX`, `PW_X`), then `$HOME/qe/bin/`, then `qe-7.5/build/bin/`, then `$PATH`. Raises a descriptive error if none resolve.
- `pwdft_validation/qe/lock.py` — `machine_lock(agent: str, desc: str)` context manager shelling out to `.claude/bin/machine-lock`. Acquires on enter, releases on exit (including on exception).
- `pwdft_validation/qe/runner.py` — `run_pwx(input_path, *, np=None, output=None, timeout=600) -> QeResult`. Auto-detects `np` via `os.cpu_count()`. Holds the lock while `mpirun -np $NP pw.x -in ...` runs.
- `pwdft_validation/qe/parse.py` — promote the regex scrapers currently duplicated in `energy.py`/`fermi.py`/`pbe_refs.py` into a single `QeScfOutput` dataclass with accessors: `total_energy_ry`, `fermi_ev`, `gamma_eigs(spin)`, `magnetization`, `n_iterations`, `converged`, `per_term_energies`, `k_weights`, `bands`. Existing diagnostics import from here.
- `pwdft_validation/qe/reference.py` — orchestrator. `regenerate(system: str | None = None)` drives `run_pwx` for every cell in a declarative `REFERENCE_CELLS` table, parses the output, and rewrites `data/qe/reference_data.toml` atomically (write to temp file, fsync, rename). Compares the fresh values against the committed TOML and exits non-zero on drift — this is the regression gate.

New CLI surface (extends the existing cyclopts app):

```text
pwdft-validate qe run <system> [--functional lda|pbe] [--np N] [--timeout 600]
pwdft-validate qe regenerate [--system si_diamond ...]   # re-writes reference_data.toml
pwdft-validate qe check                                  # parse cached .out files; no run
pwdft-validate qe compare                                # run every cell, diff vs TOML, rc!=0 on drift
```

The existing `energy components / trace / join` commands delegate to `qe.parse.QeScfOutput` internally — duplicated regexes removed.

### Phase B — Lock orchestration moves behind the Python wrapper

Every QE invocation in the repo goes through `pwdft-validate qe run ...`, which holds the machine lock internally. CLAUDE.md § Machine Coordination table gains a row:

> **QE runs:** prefer `pwdft-validate qe run <system>`. The Python harness acquires the machine lock itself and records the owning worktree — you do not need to wrap it in `.claude/bin/machine-lock run`. Nested-acquire is safe (the wrapper detects an already-held lock by the same worktree and no-ops).

The `machine-lock run cmd...` invocation stays supported for `cargo` and profiling wrappers; QE just gains a higher-level entry.

Unit test: `tests/test_qe_lock.py` spawns two concurrent `run_pwx` calls (with a no-op binary) and asserts they serialize via the shared lock directory.

### Phase C — Skill rename + rewrite: `qe-runner` → `qe-validation`

1. `git mv .claude/skills/qe-runner .claude/skills/qe-validation`.
2. Rewrite `SKILL.md`:
   - Drop the 50-line "here is the `mpirun` recipe" section.
   - Replace with "use `pwdft-validate qe run`". Include one copy-pastable example and one flowchart that routes {generate reference, diff reference, diagnose residual, run band/DOS} → concrete CLI invocation.
   - Keep the `references/pw-input.md` etc. as-is (they document input-deck format — still useful for authoring new decks).
3. Update the skill's `description:` frontmatter trigger list: add "regenerate reference", "QE regression", "validate against QE" alongside the existing triggers.
4. Global grep: `qe-runner` → `qe-validation` across `.claude/agents/*.md`, `CLAUDE.md`, `proposals/`. Redirect stale references to the new path.

### Phase D — CI job: `qe-validation` gate + per-PR accuracy/runtime report

#### D.1 — Report generator in the Python harness

Extend `pwdft_validation.qe.reference` with a report emitter:

```text
pwdft-validate qe report [--format markdown|json] [--baseline <ref.json>]
                         [--output <path>] [--fail-on-regression]
```

Behavior:

- Runs every cell in `REFERENCE_CELLS` (or the subset matching `--cells`), times the QE wall-clock per cell, then runs the pwdft-core SCF for the same cell (via a fixture YAML shipped alongside each QE `.in`), and times that.
- Parses both outputs into a `QeScfOutput` + pwdft `ScfResult` and computes per-cell accuracy deltas: `ΔE_total`, `ΔE_F`, `|Σ w_k f·ε − QE_1e|` (BSUM), magnetization delta for spin-polarized cells.
- Emits a Markdown table (the canonical on-PR format) with columns:

```text
| System | Functional | ΔE_total (meV) | ΔE_F (meV) | BSUM (meV) | QE wall | pwdft wall | Status |
|--------|-----------|---------------:|-----------:|-----------:|--------:|-----------:|:------:|
| Si     | LDA       |          12.4  |        6.4 |        0.3 |   34 s  |      11 s  |  ✅    |
| Fe BCC | LDA       |        1970.1  |      150.2 |     1200.1 |  118 s  |      52 s  |  ⚠️    |
| ...    |           |                |            |            |         |            |        |
```

- A `json` mode writes the same data as machine-readable records (one per cell) for downstream tooling (dashboards, historical trend plots).
- `--baseline <ref.json>` compares against a prior report and adds a Δ-vs-baseline column so reviewers can see whether *this PR* improved or regressed accuracy — the primary reason this is in CI.

Tolerances live in a single `TOLERANCE_TABLE` dict keyed by `(system, functional)`, shared between `qe compare` and `qe report` and mirrored by the Rust `qe_validation.rs` thresholds. Status column emoji: ✅ below tolerance; ⚠️ over tolerance but within a soft band (no fail); ❌ hard-fail tolerance (PR blocked).

#### D.2 — GitHub Actions job + PR comment

Add a `qe-validation` job to `.github/workflows/ci.yml`:

```yaml
  qe-validation:
    name: QE regression gate
    runs-on: ubuntu-latest
    if: ${{ !contains(github.event.pull_request.labels.*.name, 'docs-only') }}
    permissions:
      contents: read
      pull-requests: write   # needed to post the PR comment
    steps:
      - uses: actions/checkout@v6
      - uses: actions/cache@v4
        with:
          path: ~/conda_pkgs
          key: qe-${{ runner.os }}-7.5
      - name: Install QE 7.5 (conda-forge)
        run: |
          curl -sL https://micro.mamba.pm/install.sh | bash
          micromamba install -y -p $HOME/qe -c conda-forge quantum-espresso=7.5
          echo "$HOME/qe/bin" >> $GITHUB_PATH
      - uses: actions/cache@v4
        with:
          path: ~/.cache/pwdft-qe-baseline
          key: qe-baseline-${{ hashFiles('data/qe/reference_data.toml') }}
      - name: Install uv
        uses: astral-sh/setup-uv@v7
        with: { enable-cache: true }
      - name: Sync dependencies
        run: uv sync
      - name: Build pwdft-core (release)
        run: cargo build --release -p pwdft-core
      - name: Run QE validation + generate report
        id: report
        run: |
          uv run pwdft-validate qe report \
            --format markdown \
            --output qe-report.md \
            --baseline ~/.cache/pwdft-qe-baseline/main.json \
            --fail-on-regression
      - name: Upload JSON artifact
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: qe-report-json
          path: qe-report.json
      - name: Post / update PR comment
        if: github.event_name == 'pull_request' && always()
        uses: marocchino/sticky-pull-request-comment@v2
        with:
          header: pyqe-validation-report
          path: qe-report.md
```

Exit status:

- Report step fails (`--fail-on-regression`) if any ❌ cell exists OR if a previously-✅ cell moved to ⚠️ relative to the baseline.
- The comment is posted regardless (via `if: always()`), so reviewers see the delta table even when the job fails — the red cell is highlighted in the table.

The `sticky-pull-request-comment` action updates a single comment (`header: pyqe-validation-report`) across pushes rather than appending a new one per commit, keeping PR threads readable.

#### D.3 — Scope tiers

- **Fast tier (default on every PR)** — Si / C / Al × LDA (3 cells). Budget: ~90 s QE + ~15 s pwdft wall = under 2 min CI. Catches regressions on the structurally-green cells.
- **Full tier (nightly + on-demand via workflow_dispatch)** — all 16 cells (8 systems × LDA/PBE). Budget: ~12–15 min. Runs from a companion `.github/workflows/qe-full.yml` on `schedule: '0 3 * * *'`; posts its report to a pinned tracking issue and fails the issue comment on regression rather than blocking a PR.
- **Opt-in full on PR** — add label `qe-full` to run the full matrix on a specific PR (checked by the job's `if:` guard).

The baseline JSON cached at `~/.cache/pwdft-qe-baseline/main.json` is refreshed by a separate `main`-branch workflow after each merge (`.github/workflows/qe-baseline.yml`) so PRs diff against the current tip of main, not a drifting static snapshot.

### Phase E — Retire hand-coded invocations

Sweep + replace:

- `.claude/skills/qe-validation/SKILL.md` — done in Phase C.
- `CLAUDE.md` § Machine Coordination example — point at `pwdft-validate qe run`.
- `data/qe/README.md` — regeneration note updated to `uv run pwdft-validate qe regenerate`.
- Any proposal body referencing `mpirun pw.x -in ...` stays as historical documentation; do not rewrite completed-proposal archives.

### Ordering & dependencies

| Phase | Depends on | Unlocks |
|---|---|---|
| A | (nothing) | B, C, D |
| B | A | CI trust in the Python lock |
| C | A | future agent self-service |
| D | A, B (for lock safety in CI) | VQEF regression gate |
| E | A, C | full cleanup |

Phases A+B land together as one PR (runner + lock are a single mental unit). Phase C is a second PR (skill rename). Phase D is a third PR (CI job + workflow file). Phase E is a small cleanup PR. Parallel work possible once A lands.
