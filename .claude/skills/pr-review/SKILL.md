---
name: pr-review
description: Engineering Manager PR review helper. Fetches the PR metadata and diff and prints the EM merge checklist. Trigger on "/pr-review <PR#>".
user_invocable: true
---

# /pr-review — EM PR inspection helper

Fetch PR details for `$ARGUMENTS` (a PR number) and print the review checklist so the EM can work through it.

## Steps

1. **Parse** — `$ARGUMENTS` must be a positive integer PR number. Reject otherwise.

2. **Fetch metadata:**

   ```bash
   gh pr view "$ARGUMENTS" --json number,state,mergeable,title,author,headRefName,additions,deletions,changedFiles,body,labels
   ```

3. **Fetch the diff** (for scanning):

   ```bash
   gh pr diff "$ARGUMENTS"
   ```

4. **Present the EM checklist** (from `.claude/agents/engineering-manager.md` § PR review checklist). Tick the boxes you can verify from metadata + diff alone; flag the rest for the EM:

   - [ ] Title format: `<PROPOSAL-ID>: <description>`
   - [ ] Branch: `<PROPOSAL-ID>/<slug>`; rebased on `origin/main`
   - [ ] Changes match the proposal's Implementation section — no scope creep
   - [ ] `cargo test` passes (check CI status)
   - [ ] `cargo clippy -q --all-targets` clean
   - [ ] `cargo clippy -q --all-targets --features gpu` clean
   - [ ] `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` clean
   - [ ] Tier-2 outcome reported if the diff touches the Tier-2 trigger list
   - [ ] No new `unwrap()` or `panic!()` in production paths
   - [ ] Physics changes have verification (QE comparison, numerical test)
   - [ ] PR body has Summary + Test Plan
   - [ ] Commits follow `<ID>: <description>` format
   - [ ] **Session logbook file present** — scan the diff for `.claude/logbooks/<role>/YYYY-MM-DD-*.md`. Missing = REQUEST-CHANGES.

5. **Recommend** whether to spawn the Code Reviewer and / or Researcher:

   - **Code Reviewer** when touched LOC ≥ 200, hot-path SCF / physics, or new public-API surface. Skip for proposal-only, INDEX admin, logbook, pure-move refactors, docstring-only.
   - **Researcher** whenever the PR claims a physics bugfix or QE-validation change (run in parallel with Code Reviewer).

6. **Tier-2 trigger check** — scan the diff for paths in `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or numerics-dep bumps in `Cargo.toml`. If any hit, remind the EM to require a Tier-2 run.

7. **Report** the filled-in checklist and recommendations back to the user.

## Not in scope

- Merging the PR (that's `/merge <PR#>`).
- Opening a PR (that's `/pr-submit` for sub-agents).
- Running the quality gate yourself (that's the sub-agent's job before pushing).
