# Worktree isolation protocol

All code-editing sub-agents declare `isolation: worktree` in their frontmatter, so the harness spawns each session in a dedicated git worktree under `.claude/worktrees/agent-*`. As defense-in-depth, the `check-worktree.sh` PreToolUse hook also blocks Edit/Write/MultiEdit on paths outside the worktree (main checkout, other agents' worktrees, anything except `/tmp/`). If the hook denies a write, your `file_path` is wrong — fix the path, don't disable the hook.

The Engineering Manager is the exception: it runs on the main checkout (no `isolation` frontmatter) because it edits `proposals/` and `.claude/logbooks/engineering-manager/` directly without going through a PR.
