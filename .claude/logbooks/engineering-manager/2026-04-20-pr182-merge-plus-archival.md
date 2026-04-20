# 2026-04-20 — PR #182 merge + VGCH-2D/2E/2F archival

## Merge

**PR #182** — `VGCH followups: PZPW + CNLC + VNLM-CUD proposal drafts` — squash-merged as `5e48558`. Remote branch auto-deleted (repo has `deleteBranchOnMerge: true`). Researcher's worktree `agent-a556b5ad` removed, local branch deleted. Main-checkout rebase onto `origin/main` required because local had 2 direct-to-main commits ahead (`dab9c86 Updated agent configs`, `236c747 CI check renames`); stashed user's WIP agent-def edits around the rebase, popped cleanly.

## Archival (GRM, commit `3d6d173`)

Moved three completed VGCH follow-ups from active `proposals/` to `proposals/completed/` with `status: completed` frontmatter flipped:

- VGCH-2D — Fe LDA Class B diagnostic (PR #177)
- VGCH-2E — C diamond Class C transplant (PR #178)
- VGCH-2F — Part C session-2 ρ_core pins + V_NL cross-check (PR #179)

Flagged by the Researcher during PR #182 drafting. The three bodies had been living in active `proposals/` since their PRs merged earlier today — archival was a grooming miss the Researcher caught while reading them for physics context.

## VGCH-2F docs-drift fix (same commit)

Three references to `scripts/validate/rho_core_g_reference.py` in VGCH-2F's body were updated to point at the current `pwdft/pwdft-validation/pwdft_validation/reference/nlcc.py` location. Script migrated from `scripts/validate/` into the validation Python package post-#179 merge; VGCH-2F's body still carried the old paths. Kept the originals as "at landing time" markers so the git-archaeology cross-references stay legible.

## FLUPs carried forward (not actioned this turn)

- **Si LDA tolerance post-PZPW promotion** — if PZPW-F flips the LDA default to PW92, Si LDA may shift 30–60 meV (r_s sign-transition regime). Current Si LDA VQEF tolerance is GREEN at 0.023 meV post-SiEF-B1; may need re-baseline. Track on PZPW-F follow-up checklist.
- **Multi-functional LDA grid-lifter surface** — after PZPW Phase 1, pwdft-rs will have both PZ and PW92 LDA grid lifters. If PZPW-F promotes PW92 to default, retire the PZ grid lifters (keep `perdew_zunger_correlation` as a library helper). Not in PZPW's scope.
- **CNLC `disable_nlcc` diagnostic flag lifecycle** — remove after CNLC closes. Track on CNLC closure checklist.
- **Stashed unrelated diff in removed worktree `agent-a556b5ad`** — Researcher reported an unstaged modification to `pwdft/pwdft-core/src/symmetry/density/g_space.rs` (reorder `n_ops ≤ 1` early-return before dims asserts) when that worktree booted. Stashed locally to keep the PR clean. Worktree now removed; the stash went with it. If the change matters, it must be re-derived from git reflog on the main checkout or flagged back to whoever originated it.

## Main-checkout state after this turn

- `origin/main` at `5e48558` (post-#182).
- Local `main` at `3d6d173` (this archival commit) — ahead of origin by 3: `3d6d173 GRM archive` + `dab9c86 Updated agent configs` + `236c747 CI check renames`. Two of those are the user's direct-to-main work; the archival is EM housekeeping. User to push when ready; EM did not push.
- Uncommitted in main working tree: user's WIP agent-def edits across all six role files.

## Concurrent agents still in flight

- ROTI (Core Engineer, `agent-a7fe9ca5`) — i8→i32 rotation revert.
- TDBG (Core Engineer, `agent-a692f78b`) — CI Tier-1 `--profile=dev`.
- ESPL (Core Engineer, `agent-af025721`) — `ElectronSettings` split.

Task list tracks each one; will mark in_progress when their PRs land.

## Stale EM worktree

`agent-a45768f9` (this session's EM worktree) holds stale INDEX.md + no shipping-log.md because my earlier Write calls targeted the main-checkout absolute paths, not the worktree copies. It's on branch `worktree-agent-a45768f9` at the Researcher's commit `2d83519` (pre-squash). Leaving it for now — doesn't block anything and I'm still operating out of it via `git -C <main>` for main-affecting commands.
