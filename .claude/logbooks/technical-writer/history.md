# Technical Writer Logbook

Entries: date, what was audited, coverage stats, gaps found. Track documentation debt.

## 2026-04-16 — Handoff and orientation

Docstring conventions established: `///` with formula + parameter units in parens (Å, eV, e/ų). Gold standard: `potential/nonlocal.rs`. 13 `docs/*.md` files created replacing monolithic `theory.md` + `math-audit.md`. Index at `docs/contents.md`.

Well-documented: xc.rs, nonlocal.rs, ewald.rs, smearing.rs, pseudopotential/mod.rs.
Bare: density.rs, initial_density.rs, potentials.rs, grid.rs, context.rs, basis.rs, fft.rs, hartree.rs, local.rs, gpu/mod.rs, kpoints.rs, symmetry/*.

Open doc follow-ups from MSTR: reciprocal() docstring, Anderson constraint explanation, auto_q_tf citation, total_energy V_G0 correction. Spin exchange equivalence at xc.rs:249 needs documenting (Researcher confirmed correctness).
