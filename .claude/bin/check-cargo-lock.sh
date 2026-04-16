#!/usr/bin/env bash
# PreToolUse hook for Bash: enforce machine lock for cargo commands.
# Blocks bare cargo test/build/bench/clippy unless wrapped in machine-lock.
set -euo pipefail

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

    python3 -c "
import json, sys
print(json.dumps({
    'hookSpecificOutput': {
        'hookEventName': 'PreToolUse',
        'permissionDecision': 'deny',
        'permissionDecisionReason': sys.argv[1]
    }
}))" "Cargo commands must use the machine lock. Use: .claude/bin/machine-lock run \"<role>\" \"<desc>\" -- cargo <args>, or use the /cargo skill."
    exit 0
fi

exit 0
