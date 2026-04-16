---
name: proposal
description: Create, list, or complete engineering proposals in the proposals/ directory. Use this skill whenever the user wants to propose a change, plan an improvement, document a bug fix strategy, or archive completed work. Trigger on "/proposal", "write a proposal", "propose", "let's plan", or any mention of creating/completing proposals.
---

# Proposal Management

Proposals are the project's decision and planning records — a cross between Jira tickets and research notes. They capture not just *what* to do, but *why* it matters and the research that backs the approach. A good proposal has enough depth that someone could pick it up cold, understand the motivation, and implement it without guessing.

## ID System

Proposals use **4-letter uppercase IDs** (e.g., `SIMP`, `CBRT`, `ERRH`) instead of sequential numbers. This avoids numbering conflicts when multiple agents work concurrently. IDs should be mnemonic — easy to remember and clearly related to the topic.

Files are named `XXXX-slug.md` (e.g., `SIMP-simpson-radial-quadrature.md`).

The index at `proposals/INDEX.md` lists all active and completed proposals with their IDs.

## Subcommands

Parse the user's argument to determine which operation:

- `create <topic>` — write a new proposal (default if just a topic is given)
- `list` — show open and completed proposals
- `complete <id>` — archive a finished proposal

If the argument is ambiguous, assume `create`.

---

## `create <topic>`

### 1. Check for duplicates and overlap

**This step is mandatory.** Before writing anything:

- Read `proposals/INDEX.md` to see all active and completed proposals.
- Scan active proposal files in `proposals/` for overlapping scope — read titles and Problem sections.
- Scan completed proposals in `proposals/completed/` — has this been tried before?
- If there is overlap with an existing active proposal, **do not create a new one**. Instead, report the overlap and ask whether to extend the existing proposal or proceed anyway.
- If the topic was previously completed, flag it and ask whether it needs reopening.

### 2. Choose a 4-letter ID

Pick a mnemonic 4-letter uppercase ID. Verify it doesn't collide with any existing ID in `proposals/INDEX.md` (both active and completed sections). Good IDs are abbreviations of the core concept: `SIMP` for Simpson's rule, `CBRT` for cbrt optimization, `ERRH` for error handling.

### 3. Research before writing

This is critical. Proposals must be grounded in the actual codebase, not speculation. Before writing a single line of the proposal:

- **Read the relevant source files.** If proposing a change to mixing, read `src/scf/mixing.rs`. If proposing a new lint, run `cargo clippy` and count real warnings.
- **Quantify the problem.** Use grep/glob to count occurrences, measure actual impact, find all affected locations. Put real numbers in the proposal (e.g., "25 unwrap() calls remain" not "many unwrap() calls").
- **Validate feasibility.** If the proposal depends on a crate or API, verify it exists and works. If it claims a function signature needs changing, read the function and confirm.

Do NOT write proposals based on assumptions. Every claim should be verifiable by reading the code or running a command.

### 4. Write the proposal

Create `proposals/XXXX-slug.md` where `XXXX` is the 4-letter ID and `slug` is a short kebab-case summary.

Use this structure:

```markdown
---
id: XXXX
status: active
priority: critical|high|medium|low
complexity: trivial|small|medium|large
risk: low|medium|high
depends_on: []
blocks: []
---

# XXXX: Title

## Problem

What's wrong or missing, with evidence. Include counts, file paths, and concrete
examples from the codebase. If there's a bug, show how to reproduce it. If there's
a performance issue, show a measurement. This section answers "why should we care?"

## Research

(Optional — include when alternatives were considered or analysis was performed.)

Analysis, measurements, alternatives considered, and why the chosen approach wins.
Tables comparing options are encouraged. Link to external references (papers, docs)
where relevant.

## Implementation

Concrete steps with code examples where helpful. Each step should name specific
files and functions being changed. Use code blocks for non-trivial changes.
Order steps by dependency (what must happen first) or by risk (highest-impact first).

## Verification

How to confirm the implementation is correct. This could include:
- Specific tests to write or run
- Cargo commands to validate (clippy, test, bench)
- Quantum ESPRESSO comparisons for physics changes
- Before/after measurements for performance changes
```

#### Frontmatter field guide

| Field | Values | Guidance |
|-------|--------|----------|
| `status` | `active`, `deferred`, `completed` | Use `deferred` for future work not yet prioritized |
| `priority` | `critical`, `high`, `medium`, `low` | `critical` = blocks research correctness; `high` = foundation for other work; `medium` = valuable improvement; `low` = nice to have |
| `complexity` | `trivial`, `small`, `medium`, `large` | `trivial` = under 1hr, single file; `small` = 1-2hrs, few files; `medium` = half day; `large` = full day+ |
| `risk` | `low`, `medium`, `high` | `low` = no physics impact; `medium` = touches physics code; `high` = could break convergence or results |
| `depends_on` | list of IDs | Proposals that must be completed first |
| `blocks` | list of IDs | Proposals that cannot start until this one completes |

#### Style guidance

- Lead with data, not opinion. "15 float comparisons use `==`" beats "we should probably fix float comparisons."
- Use tables to organize findings.
- Include code blocks for any non-trivial change — show the before/after.
- Keep the tone technical but direct. These are working documents, not formal specs.
- The Research section is optional — skip it if the change is straightforward.

### 5. Update the index

Add the new proposal to `proposals/INDEX.md` in the appropriate priority section. One row with ID, title, complexity, risk, depends_on, and blocks.

### 6. Present to the user

After writing, give a brief summary: the ID, title, problem statement in one sentence, and the key implementation steps. Ask if they want to adjust anything before it's finalized.

---

## `list`

Read and display `proposals/INDEX.md`. It is organized by priority tier with all metadata visible.

If the index seems stale (files exist that aren't listed, or listed files don't exist), note the discrepancies.

---

## `complete <id>`

Archive a proposal after implementation is done and verified.

### 1. Verify tests pass

Run `cargo test` (or `cargo test --features gpu` if the proposal touches GPU code). Do NOT proceed if tests fail — report the failures and stop.

### 2. Run clippy

Run `cargo clippy -q --all-targets`. Report any new warnings introduced by the implementation.

### 3. Move the file

Move `proposals/XXXX-*.md` to `proposals/completed/XXXX-*.md`.

### 4. Update the index

Move the proposal's row from the Active section of `proposals/INDEX.md` to the Completed section. Update `status` in the file's frontmatter to `completed`.

### 5. Update dependents

Check if any other active proposals had this ID in their `depends_on` list. If so, remove it (the dependency is now satisfied). Note any proposals that are now unblocked.

### 6. Confirm

Report that the proposal has been archived, with a summary of what was implemented, that tests passed, and which proposals (if any) are now unblocked.
