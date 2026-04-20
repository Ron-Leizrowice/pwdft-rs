# 2026-04-20 — TDBG merge + no-backcompat policy landing

## Shipped

- **PR #183 TDBG merged.** Squash-SHA: see `git log origin/main --oneline`. Measured CI wall-time: test step 343 s → 246 s (−28%), full job 643 s → 315 s (−51%) on the PR's own warm-cache steady-state runs. 330 s Tier-1 budget guardrail inline in the workflow step (not the proposal's 240 s — real floor is 281 s, so 330 s keeps the pre-TDBG 343 s baseline as the regression ceiling). Proposal archived to `proposals/completed/TDBG-*.md`; INDEX row removed; shipping-log's **CI** bucket bumped to include TDBG; cumulative PR count 72 → 73.

## Policy: no backwards compatibility

User directive: "We do not endorse legacy code or backwards compatibility. This code is in active development, it has no users yet. We do not do backwards compatibility, we move forwards."

Landed in commit `ab781af`:

- New `.claude/agents/shared/no-backcompat.md` — canonical rule, forbidden-constructs list (`#[serde(alias)]`, `#[deprecated]`, "one-release warning" paths, feature-flag preservation, `if use_old_path` branches, alias re-exports), acceptance criterion (old form is a hard compile- or parse-error), and a `REQUEST-CHANGES` trigger for reviewers.
- Wired into all 6 agent defs (EM + code-reviewer + core-engineer + performance-engineer + researcher + technical-writer).
- Code-reviewer `## Mindset` section got an explicit bullet: legacy-shim PRs are an automatic `REQUEST-CHANGES` blocker.
- EM memory: `feedback_no_backcompat.md` + MEMORY.md index entry.

## ESPL proposal amendment

Landed in the same commit `ab781af`:

- **New Part C** — `ElectronsPhysics.nspin: usize` → `spin_polarized: bool`. Argued in the proposal: `nspin` was never a count; `usize` invites nonsense values (0, 3, 17); non-collinear DFT is a separate code path, not a knob. Migration table maps `nspin == 1` ↔ `!spin_polarized`; no `nspin()` accessor method.
- **Stripped legacy YAML support** from the Migration and Acceptance sections (per the no-backcompat policy). Any YAML writing `nspin:` at all is a hard parse error naming `spin_polarized`. In-repo fixtures migrate in the same PR.
- **Recovery-agent heads-up** pinned at the top of "Previous attempt" noting Part C was added after the `1ae1416` checkpoint, so the recovery agent layers it on top.

**Open question for next EM session.** The ESPL agent (`af0257216752d1873`) was already running when Part C was added and the shim was stripped. Its PR (branch `ESPL/electrons-settings-split`, WIP checkpoint PR #185) will likely return with Part A+B + a `#[serde(alias)]` shim matching the *pre-amendment* proposal. On review: either (a) merge as partial and file a follow-up for Part C + shim-strip, or (b) REQUEST-CHANGES pointing at the updated proposal and the new `shared/no-backcompat.md`. Option (b) is cleaner given the no-backcompat policy is a hard rule.

## Other in-flight work

- **ROTI** (PR #184) — WIP checkpoint only, original agent aborted. Branch `ROTI/revert-symm-rotation-i8` @ `8b6fbbc` preserved; proposal has a "Previous attempt" recovery section. Needs a fresh core-engineer to run `/quality-gate` + `/test --tier2` + `/pr-submit`. Not restarted this session.
- **Worktree cleanup** — `agent-a692f78b` (TDBG) removed as part of the merge trilogy. `agent-a7fe9ca5` (ROTI) and `agent-af025721` (ESPL, still running) remain.

## Merge trilogy note

Rebase failed the first time because of uncommitted proposal WIP (ROTI + TDBG proposal files with recovery / landed-PR sections). Solved with `git stash push -- <paths>` + `rebase origin/main` + `git stash pop`. The `/merge` skill already documents this stash/pop plumbing — followed it and it worked.

## Cruft

- Global `~/.gitignore` excludes `.claude`; new files under `.claude/` require `git add -f`. Noted for next session — the repo-local `.gitignore` has negations for specific tracked paths but new untracked files (agent memory, shared protocol docs) still need `-f`.
