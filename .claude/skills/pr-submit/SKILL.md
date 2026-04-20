---
name: pr-submit
description: Run the quality gate, rebase on origin/main, push, and open (or promote to ready) a PR against main. Use at end of implementation. Trigger on "/pr-submit".
user_invocable: true
---

# /pr-submit — Quality gate + rebase + PR

Finalize an implementation branch and either open the PR (if none exists) or promote an existing `/pr-draft` checkpoint to ready-for-review. The argument is optional; if present, it's appended to the PR body's Summary section.

Use `/pr-draft` for mid-implementation checkpoints; `/pr-submit` is the end-of-work gate.

## Steps

1. **Verify worktree + branch.** `pwd` must be under `.claude/worktrees/agent-*`; current branch must be `<PROPOSAL-ID>/<slug>`. Abort if on `main` or a stale branch.

2. **Confirm the logbook entry exists.** Check that a file `.claude/logbooks/<your-role>/YYYY-MM-DD-<slug>.md` has been created inside your worktree (see `.claude/agents/shared/session-end.md`). If missing, stop and tell the caller to write it — the entry must land with the PR, not be deferred.

3. **Run the quality gate** under the machine lock. All five checks must pass:

   ```bash
   .claude/bin/machine-lock run "<role>" "pr-submit quality gate" -- bash -c '
     cargo clippy -q --fix --allow-dirty --allow-staged --all-targets &&
     cargo clippy -q --all-targets &&
     cargo clippy -q --all-targets --features gpu &&
     RUSTDOCFLAGS="-D warnings" cargo doc --no-deps &&
     cargo test'
   ```

   Replace `<role>` with your agent role. If any check fails, stop and report — do not push.

   The `cargo clippy --fix` pass may modify files. After the gate finishes, stage those fixes along with your session logbook entry and any agent-memory writes, then commit if anything is staged:

   ```bash
   git -C "$(pwd)" add .
   if ! git -C "$(pwd)" diff --cached --quiet; then
     git -C "$(pwd)" commit -m "<ID>: quality-gate fixes + session notes"
   fi
   ```

   `git add .` from the worktree root catches everything — the worktree is isolated, so every file under it belongs to this PR. This is the point at which memory and logbook writes enter the PR; don't skip it.

4. **Tier-2 check.** If the diff touches `pwdft/pwdft-core/src/{scf,potential,symmetry,pseudopotential,eigensolver,gpu}/`, `basis.rs`, `fft.rs`, `ewald.rs`, `crystal.rs`, `kpoints.rs`, or bumps a numerics dep, also run `/test --tier2`. Capture the outcome — it must appear in the PR body's Test Plan.

5. **Rebase on `origin/main`.** Conflicts get resolved in the worktree:

   ```bash
   git -C "$(pwd)" fetch origin
   git -C "$(pwd)" rebase origin/main
   ```

6. **Push** (force-with-lease is safe after a rebase):

   ```bash
   git -C "$(pwd)" push --force-with-lease -u origin "$(git -C "$(pwd)" branch --show-current)"
   ```

7. **Read the proposal file** so you can cite it in the PR body. Extract the `<ID>` and title from `proposals/<ID>-*.md`.

8. **Open the PR or promote the existing draft.**

   ```bash
   existing=$(gh pr list --head "$(git -C "$(pwd)" branch --show-current)" --state open --json number,isDraft --jq '.[0]')
   ```

   - **If `$existing` is empty** (no `/pr-draft` checkpoint was opened earlier): `gh pr create` with the full ready body:

     ```bash
     gh pr create --title "<ID>: <description>" --body "$(cat <<'EOF'
     ## Proposal

     <ID>: <title>

     ## Summary

     - <what changed>
     - <why>

     ## Test plan

     - `cargo test` — <PASS/FAIL>
     - `cargo clippy -q --all-targets` — <PASS/FAIL>
     - `cargo clippy -q --all-targets --features gpu` — <PASS/FAIL>
     - `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps` — <PASS/FAIL>
     - Tier 2 (if touched): <N/A or outcome>
     EOF
     )"
     ```

   - **If `$existing` is a draft**: promote to ready and update the body with the final Test Plan:

     ```bash
     num=$(echo "$existing" | jq -r .number)
     gh pr edit  "$num" --title "<ID>: <description>" --body "<full ready body with Test Plan>"
     gh pr ready "$num"
     ```

   - **If `$existing` is already ready**: just update the body to reflect the final Test Plan. Do not open a duplicate.

9. **Report** the PR URL, number, branch, and the logbook filename to the user. If you promoted a draft, note that the CI runs from step 6's push should be current.

## Errors

- Quality gate fails → stop. Fix the issue, commit, and re-run `/pr-submit`. Do not skip checks.
- Rebase conflicts → resolve in the worktree, commit the resolution, re-run `/pr-submit` from step 4.
- `gh pr create` fails on "base branch identical" → the branch has no commits beyond `origin/main`; nothing to PR.

## See also

- `.claude/skills/pr-draft/SKILL.md` — mid-implementation checkpoint skill. Call this first, early and often; then `/pr-submit` at the end.
- `.claude/agents/shared/quality-gate.md`
- `.claude/agents/shared/worktree.md` § Rebase before submitting the PR
