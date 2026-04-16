#!/usr/bin/env bash
# PreToolUse hook for Edit/Write: block source file edits in the main checkout.
# Source files (src/, tests/, benches/) must be edited in worktrees.
# Infrastructure files (.claude/, proposals/, CLAUDE.md, docs/, etc.) are allowed.
set -euo pipefail

input=$(cat)
file_path=$(python3 -c "import sys,json; print(json.load(sys.stdin).get('tool_input',{}).get('file_path',''))" <<< "$input" 2>/dev/null || echo "")

[ -z "$file_path" ] && exit 0

# Determine repo root
toplevel=$(git rev-parse --show-toplevel 2>/dev/null || echo "")
[ -z "$toplevel" ] && exit 0

# In a worktree, .git is a file; in the main checkout, it's a directory
[ -f "$toplevel/.git" ] && exit 0

# We're in the main checkout. Check if this is a source file.
rel_path="${file_path#$toplevel/}"
case "$rel_path" in
    src/*|tests/*|benches/*|build.rs)
        python3 -c "
import json, sys
print(json.dumps({
    'hookSpecificOutput': {
        'hookEventName': 'PreToolUse',
        'permissionDecision': 'deny',
        'permissionDecisionReason': sys.argv[1]
    }
}))" "Source file '${rel_path}' must be edited in a worktree, not the main checkout. Use EnterWorktree (interactive) or isolation: \"worktree\" (Agent tool)."
        exit 0
        ;;
esac

exit 0
