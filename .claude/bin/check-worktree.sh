#!/usr/bin/env bash
# PreToolUse hook for Edit/Write/MultiEdit: enforce write isolation.
#
# Two rules, applied in order:
#
#   1. WORKTREE WRITE ISOLATION
#      If the *current shell* is in a git worktree (i.e. the cwd's toplevel
#      has .git as a *file*, which is how worktrees mark themselves), then
#      every Edit/Write target MUST resolve inside that same worktree.
#
#      This prevents the most common isolation leak: a sub-agent in
#      /Users/.../.claude/worktrees/agent-X passing an absolute path like
#      /Users/.../src/foo.rs (which lands in the main checkout) to Edit.
#      Without this rule the previous version of the hook checked the cwd's
#      worktree status and approved the write because cwd was inside a
#      worktree, even though the *target* was outside it.
#
#   2. MAIN-CHECKOUT SOURCE-FILE GUARD
#      If the cwd is the main checkout (toplevel's .git is a *directory*),
#      block edits to src/, tests/, benches/, and build.rs. Allow
#      .claude/, proposals/, CLAUDE.md, docs/, etc. — those are
#      infrastructure that may legitimately be changed from main as admin
#      commits.
#
# /tmp is always allowed (scratch / rescue).
#
# Notes:
# - We use python3 for path realpath comparison because bash's `realpath`
#   isn't on macOS by default.
# - We swallow git failures silently (allow) — if we're not in any repo,
#   nothing to enforce here.
set -euo pipefail

input=$(cat)
file_path=$(python3 -c "import sys,json; print(json.load(sys.stdin).get('tool_input',{}).get('file_path',''))" <<< "$input" 2>/dev/null || echo "")

[ -z "$file_path" ] && exit 0

# Always allow /tmp scratch.
case "$file_path" in
    /tmp/*) exit 0 ;;
esac

# Determine the cwd's repo toplevel. If we're not in a repo, allow.
cwd_toplevel=$(git rev-parse --show-toplevel 2>/dev/null || echo "")
[ -z "$cwd_toplevel" ] && exit 0

deny() {
    python3 -c "
import json, sys
print(json.dumps({
    'hookSpecificOutput': {
        'hookEventName': 'PreToolUse',
        'permissionDecision': 'deny',
        'permissionDecisionReason': sys.argv[1]
    }
}))" "$1"
    exit 0
}

# Realpath helper using python (portable across macOS/Linux).
realpath_py() {
    python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$1"
}

if [ -f "$cwd_toplevel/.git" ]; then
    # === Rule 1: shell is inside a worktree ===
    # The file path MUST resolve to somewhere under this worktree's root.
    abs_file=$(realpath_py "$file_path")
    abs_root=$(realpath_py "$cwd_toplevel")
    case "$abs_file" in
        "$abs_root"|"$abs_root"/*)
            exit 0
            ;;
        *)
            deny "Worktree isolation violation: target '${file_path}' resolves to '${abs_file}', which is outside the current worktree '${abs_root}'. Sub-agents must only Edit/Write/MultiEdit files inside their assigned worktree. If you need to reference data from elsewhere, copy it into your worktree or use /tmp scratch."
            ;;
    esac
fi

# === Rule 2: shell is in the main checkout ===
# Block source-file edits (must be done in a worktree via PR), allow infra.
rel_path="${file_path#$cwd_toplevel/}"
case "$rel_path" in
    src/*|tests/*|benches/*|build.rs)
        deny "Source file '${rel_path}' must be edited in a worktree, not the main checkout. Use EnterWorktree (interactive) or isolation: \"worktree\" (Agent tool)."
        ;;
esac

exit 0
