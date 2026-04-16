---
name: cargo
description: Run cargo commands with automatic machine lock coordination. Use /cargo instead of bare cargo to prevent interfering with other agents' benchmarks or tests. Trigger on "/cargo".
user_invocable: true
---

# /cargo — Machine-coordinated cargo command

Run the cargo command `$ARGUMENTS` with automatic machine lock acquire/release.

## Steps

1. Determine your agent role. If you are a specific role (Core Engineer, Performance Engineer, Researcher, Code Reviewer, Technical Writer), use that. Otherwise use "Agent".

2. Run this single command (it acquires the lock, runs cargo, and releases automatically):

```bash
.claude/bin/machine-lock run "<your-role>" "cargo $ARGUMENTS" -- cargo $ARGUMENTS
```

3. If the lock is **BLOCKED**, report which agent holds it and what they're doing. Do NOT force-release. Suggest waiting and retrying.

4. Report the cargo output to the user. If the command failed, report the error.

## Examples

- `/cargo test` → `.claude/bin/machine-lock run "Core Engineer" "cargo test" -- cargo test`
- `/cargo bench --bench scf_benchmarks` → `.claude/bin/machine-lock run "Performance Engineer" "cargo bench" -- cargo bench --bench scf_benchmarks`
- `/cargo clippy -q --all-targets` → `.claude/bin/machine-lock run "Code Reviewer" "cargo clippy" -- cargo clippy -q --all-targets`
