# Worktree isolation protocol

All sub-agent sessions run inside a dedicated git worktree under `.claude/worktrees/agent-*`. The `check-worktree.sh` PreToolUse hook blocks Edit/Write/MultiEdit on paths outside the worktree (main checkout, other agents' worktrees, anything else except `/tmp/`). If the hook denies a write, your `file_path` is wrong — fix the path, don't disable the hook.

## At session start

```bash
pwd                 # MUST resolve to .claude/worktrees/agent-*
git worktree list   # confirm the branch is checked out where you think
```

If `pwd` is the main checkout, stop and report a harness failure.

## Branch from `origin/main`, not local `main`

```bash
git -C "$(pwd)" fetch origin
git -C "$(pwd)" checkout -b <PROPOSAL-ID>/<slug> origin/main
```

Local `main` may lag; `origin/main` is the source of truth.

## Ground rules

- Use `git -C "$(pwd)"` for every git command so Bash cwd drift can't pin operations to the wrong tree.
- Never write absolute paths starting at the main checkout root — those resolve to the main checkout and the hook will block. Use relative paths, or paths that begin with your worktree root.
- Read access to the main checkout is fine (proposals, source, CLAUDE.md). Write access is not.
- `/tmp/` is always writable for scratch.

## Rebase before submitting the PR

```bash
git -C "$(pwd)" fetch origin
git -C "$(pwd)" rebase origin/main
git -C "$(pwd)" push --force-with-lease origin <branch>
```

Resolve conflicts in your worktree. This prevents "DIRTY/CONFLICTING" PRs that force the EM to rebase manually.
