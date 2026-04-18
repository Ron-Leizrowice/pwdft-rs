#!/usr/bin/env bash
# PreToolUse hook for Bash: enforce machine lock + owner scope for cargo/QE.
#
# Rules (evaluated in order):
#   1. Command doesn't look like a cargo / heavy CPU job → allow.
#   2. Command is self-wrapping in `machine-lock` → allow (the script acquires
#      atomically before running).
#   3. A new-format lock is held:
#        - if the current Bash call's cwd is inside the lock's recorded
#          worktree root, allow (the owning agent running further commands);
#        - otherwise deny with an informative message.
#   4. A legacy-format lock exists → treat as stale, deny. Ask the caller to
#      wrap via `machine-lock run` (which will clear it).
#   5. No lock → deny and instruct the caller to use machine-lock.
#
# MLFX fix (2026-04-18): rule 3 is new. Pre-MLFX the hook allowed *any*
# cargo command whenever the lockfile existed, so one agent holding the
# lock during a bench could be interrupted by another agent running
# `cargo test` and saturating the CPU (observed in WFRX #99).
set -euo pipefail

# Resolve the main repo (shared across worktrees) so we look at the same
# lock directory the `machine-lock` script writes to.
COMMON_GIT_DIR="$(git rev-parse --git-common-dir 2>/dev/null || echo "")"
if [ -n "$COMMON_GIT_DIR" ]; then
    COMMON_GIT_DIR="$(cd "$COMMON_GIT_DIR" && pwd)"
    MAIN_REPO="$(dirname "$COMMON_GIT_DIR")"
else
    MAIN_REPO="$(cd "$(dirname "$0")/../.." && pwd)"
fi
LOCKBASE="$MAIN_REPO/.claude/locks/machine.lock"
LOCKD="$LOCKBASE.d"

input=$(cat)

# Extract command and cwd from the PreToolUse JSON input.
command=$(python3 -c "import sys,json; print(json.load(sys.stdin).get('tool_input',{}).get('command',''))" <<< "$input" 2>/dev/null || echo "")
tool_cwd=$(python3 -c "
import sys, json, os
doc = json.load(sys.stdin)
cwd = doc.get('tool_input', {}).get('cwd') or doc.get('cwd') or os.environ.get('PWD', '')
print(cwd)
" <<< "$input" 2>/dev/null || echo "")
[ -z "$tool_cwd" ] && tool_cwd="${PWD:-$(pwd)}"

[ -z "$command" ] && exit 0

# Match cargo compile/test/bench/clippy/run at command position.
# Match: "cargo test", "cd foo && cargo test", "X ; cargo bench"
# Don't match: "echo 'cargo test'", "grep cargo", etc.
if ! echo "$command" | grep -qE '(^|&&|\|\||;)\s*cargo\s+(test|build|bench|clippy|run)'; then
    exit 0
fi

# Allow if the command wraps itself in machine-lock (acquire will happen
# atomically, including --wait inside `run`).
if echo "$command" | grep -qF 'machine-lock'; then
    exit 0
fi

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

_realpath() {
    python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$1"
}

# Rule 4: legacy flat-file lock → stale, deny. Don't auto-clear here — we're
# a hook, not a manager. Tell the user to run through machine-lock.
if [ -f "$LOCKBASE" ] && [ ! -d "$LOCKD" ]; then
    deny "Legacy-format lock detected at ${LOCKBASE}. Run '.claude/bin/machine-lock release' to clear it, then '.claude/bin/machine-lock run ...' to run your cargo command."
fi

# Rule 3: new-format lock exists — check ownership by worktree root.
if [ -d "$LOCKD" ]; then
    lock_agent=$(cat "$LOCKD/agent" 2>/dev/null || echo "unknown")
    lock_desc=$(cat "$LOCKD/desc" 2>/dev/null || echo "unspecified")
    lock_worktree=$(cat "$LOCKD/worktree" 2>/dev/null || echo "")
    lock_pid=$(cat "$LOCKD/pid" 2>/dev/null || echo "0")

    # If metadata is missing, treat as corrupt → deny and ask the caller to
    # release. Avoid auto-clearing here; that's the acquirer's responsibility.
    if [ -z "$lock_worktree" ]; then
        deny "Corrupt lock at ${LOCKD} (missing worktree metadata). Run '.claude/bin/machine-lock release' to clear, then retry with 'machine-lock run'."
    fi

    # PID liveness: if the owner is dead, treat the lock as stale for scoping
    # purposes. We still don't auto-clear — deny and point at machine-lock,
    # which will clear on next `acquire`.
    if [ -n "$lock_pid" ] && [ "$lock_pid" != "0" ] && ! kill -0 "$lock_pid" 2>/dev/null; then
        deny "Stale lock at ${LOCKD} (owner pid ${lock_pid} dead). Run '.claude/bin/machine-lock run \"<role>\" \"<desc>\" -- <cmd>' — it will auto-clear on acquire."
    fi

    abs_cwd="$(_realpath "$tool_cwd" 2>/dev/null || echo "$tool_cwd")"
    abs_root="$(_realpath "$lock_worktree" 2>/dev/null || echo "$lock_worktree")"
    case "$abs_cwd" in
        "$abs_root"|"$abs_root"/*)
            # Same worktree — owner is running further commands. Allow.
            exit 0
            ;;
        *)
            deny "Cargo command denied: machine lock is held by '${lock_agent}' in worktree '${abs_root}' (task: ${lock_desc}). Your cwd is '${abs_cwd}', which is outside that worktree. Wait for release, or run from your own worktree after acquiring the lock."
            ;;
    esac
fi

# Rule 5: no lock at all — tell the caller to use machine-lock.
deny "Cargo commands must use the machine lock. Use: .claude/bin/machine-lock run \"<role>\" \"<desc>\" -- cargo <args>, or acquire the lock first with: .claude/bin/machine-lock acquire \"<role>\" \"<desc>\""
