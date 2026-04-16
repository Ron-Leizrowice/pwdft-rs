#!/usr/bin/env bash
# PreToolUse hook for Bash: enforce machine lock for cargo commands.
# Blocks bare cargo test/build/bench/clippy unless:
#   1. The command is wrapped in machine-lock run, OR
#   2. The lock is already held (two-step acquire/run/release workflow)
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo "$(cd "$(dirname "$0")/../.." && pwd)")"
LOCKFILE="$REPO_ROOT/.claude/locks/machine.lock"

input=$(cat)
command=$(python3 -c "import sys,json; print(json.load(sys.stdin).get('tool_input',{}).get('command',''))" <<< "$input" 2>/dev/null || echo "")

[ -z "$command" ] && exit 0

# Check if this involves a cargo compile/test command at command position
# Match: "cargo test", "cd foo && cargo test", "X ; cargo bench"
# Don't match: "echo 'cargo test'", "grep cargo", etc.
if echo "$command" | grep -qE '(^|&&|\|\||;)\s*cargo\s+(test|build|bench|clippy|run)'; then
    # Allow if wrapped in machine-lock
    if echo "$command" | grep -qF 'machine-lock'; then
        exit 0
    fi

    # Allow if lock is already held (two-step acquire/run/release workflow)
    if [ -f "$LOCKFILE" ]; then
        exit 0
    fi

    python3 -c "
import json, sys
print(json.dumps({
    'hookSpecificOutput': {
        'hookEventName': 'PreToolUse',
        'permissionDecision': 'deny',
        'permissionDecisionReason': sys.argv[1]
    }
}))" "Cargo commands must use the machine lock. Use: .claude/bin/machine-lock run \"<role>\" \"<desc>\" -- cargo <args>, or acquire the lock first with: .claude/bin/machine-lock acquire \"<role>\" \"<desc>\""
    exit 0
fi

exit 0
