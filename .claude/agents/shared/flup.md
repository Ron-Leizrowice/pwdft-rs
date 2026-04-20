# Flagged-for-follow-up protocol

If during your session you notice work outside your role's competency, do **not** solve it and do **not** expand your PR's scope. The right action is to flag it back to the Engineering Manager so the right specialist picks it up.

## In your final return summary

Add a `## Flagged for follow-up` section listing each finding, one line per item:

```text
## Flagged for follow-up
- pwdft/pwdft-core/src/scf/mod.rs:430 — needless allocation in spin loop; Performance Engineer.
- pwdft/pwdft-core/src/potential/local.rs:78 — manual `for ... { x += ... }` sum should be `iter().sum()`; Code Reviewer.
- pwdft/pwdft-core/src/potential/xc.rs:60 — docstring says "Hartree" but function returns Rydberg; Researcher should confirm.
```

Each line: `file:line — one-sentence finding; target role.`

## Role routing

- Physics / math correctness question → **Researcher**
- Hot-path inefficiency, profile-worthy → **Performance Engineer**
- New feature, bug fix, refactor → **Core Engineer**
- Dead code, naming, idioms, lint → **Code Reviewer**
- Missing / stale docstring, README, comment → **Technical Writer**

## What the EM does with these

The EM triages each flagged item against the existing backlog. Genuinely new findings become stub proposals (4/5-letter ID, 1-paragraph problem statement) in `proposals/`. Duplicates get appended to the existing proposal's notes. See the Engineering Manager role for the FLUP seeding pattern and the "promote on signal" rule.
