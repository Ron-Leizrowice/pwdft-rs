#!/usr/bin/env bash
# machine-lock test suite (MLFX)
#
# Exercises the four MLFX bug fixes:
#   1. Hook owner check — cross-worktree cargo is denied.
#   2. Atomic acquire — two racing `machine-lock acquire` calls can't both win.
#   3. --wait mode — `acquire --wait` blocks, succeeds after release; `run`
#      serializes concurrent calls.
#   4. PID liveness staleness — a dead PID clears the lock; a live PID does
#      not. Legacy flat-file locks are treated as stale.
#
# Usage: bash .claude/bin/tests/machine-lock.test.sh
# Exits non-zero on any test failure.
#
# All state is kept in a private fixture dir under /tmp so the real
# .claude/locks/ is untouched.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LOCK="$SCRIPT_DIR/machine-lock"
HOOK="$SCRIPT_DIR/check-cargo-lock.sh"

FIXTURE_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/mlfx-test.XXXXXX")"
trap 'rm -rf "$FIXTURE_ROOT" /tmp/mlfx-*.out 2>/dev/null || true' EXIT

PASS=0
FAIL=0
FAIL_NAMES=()

pass() { PASS=$((PASS + 1)); printf '  PASS  %s\n' "$1"; }
fail() { FAIL=$((FAIL + 1)); FAIL_NAMES+=("$1"); printf '  FAIL  %s -- %s\n' "$1" "$2"; }

# --------------------------------------------------------------------------
# Fixture: build a fake main repo + two worktrees, all with their own .git
# so git rev-parse resolves the same way it would in the real agent setup.
# --------------------------------------------------------------------------

setup_fixture() {
    MAIN="$FIXTURE_ROOT/main"
    mkdir -p "$MAIN/.claude/bin" "$MAIN/.claude/locks"
    (cd "$MAIN" && git init -q && git -c user.email=t@t -c user.name=t commit --allow-empty -q -m init) 2>/dev/null

    WT_A="$FIXTURE_ROOT/wt-a"
    (cd "$MAIN" && git worktree add -q -b test-branch-a "$WT_A") 2>/dev/null

    WT_B="$FIXTURE_ROOT/wt-b"
    (cd "$MAIN" && git worktree add -q -b test-branch-b "$WT_B") 2>/dev/null

    cp "$LOCK" "$MAIN/.claude/bin/machine-lock"
    cp "$HOOK" "$MAIN/.claude/bin/check-cargo-lock.sh"
    chmod +x "$MAIN/.claude/bin/machine-lock" "$MAIN/.claude/bin/check-cargo-lock.sh"

    export TEST_LOCK="$MAIN/.claude/bin/machine-lock"
    export TEST_HOOK="$MAIN/.claude/bin/check-cargo-lock.sh"
    export TEST_LOCKDIR="$MAIN/.claude/locks"
}

reset_locks() {
    rm -rf "$TEST_LOCKDIR/machine.lock" "$TEST_LOCKDIR/machine.lock.d"
}

# Build a PreToolUse JSON payload using a Python heredoc to avoid bash-brace
# headaches. Echoes the JSON to stdout.
build_hook_input() {
    local cmd="$1" cwd="$2"
    CMD="$cmd" CWD="$cwd" python3 <<'PY'
import json, os
print(json.dumps({
    "tool_input": {"command": os.environ["CMD"], "cwd": os.environ["CWD"]},
    "cwd": os.environ["CWD"],
}))
PY
}

# Feed synthetic JSON to the hook in a cwd; return "allow" on empty output or
# "deny:<reason>" on denial.
hook_decision() {
    local cwd="$1" cmd="$2"
    local input out reason
    input=$(build_hook_input "$cmd" "$cwd")
    out=$(cd "$cwd" && printf '%s' "$input" | bash "$TEST_HOOK" 2>&1)
    if [ -z "$out" ]; then
        echo "allow"
        return
    fi
    reason=$(OUT="$out" python3 <<'PY'
import json, os, sys
try:
    doc = json.loads(os.environ["OUT"])
    print(doc.get("hookSpecificOutput", {}).get("permissionDecisionReason", "(no reason)"))
except Exception as e:
    print(f"parse-error: {e}")
PY
)
    echo "deny:$reason"
}

# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

setup_fixture

echo "==> Bug 1: hook owner check (cross-worktree denial)"

# Use a background sleep as an explicit "live owner PID" so the hook's
# liveness check is deterministic — not dependent on whatever PPID the
# test harness happens to have.
sleep 600 &
LIVE_PID=$!
# shellcheck disable=SC2329  # invoked via EXIT trap below
cleanup() {
    kill "$LIVE_PID" 2>/dev/null || true
    wait "$LIVE_PID" 2>/dev/null || true
    rm -rf "$FIXTURE_ROOT" /tmp/mlfx-*.out 2>/dev/null || true
}
trap cleanup EXIT

reset_locks
(cd "$WT_A" && ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire "Agent A" "scoped test" >/dev/null)

dec=$(hook_decision "$WT_A" "cargo test")
case "$dec" in
    allow) pass "hook allows cargo from owning worktree" ;;
    *)     fail "hook allows cargo from owning worktree" "got: $dec" ;;
esac

dec=$(hook_decision "$WT_B" "cargo test")
case "$dec" in
    deny:*worktree*) pass "hook denies cross-worktree cargo" ;;
    *)               fail "hook denies cross-worktree cargo" "got: $dec" ;;
esac

dec=$(hook_decision "$MAIN" "cargo test")
case "$dec" in
    deny:*worktree*|deny:*outside*) pass "hook denies cargo from main checkout" ;;
    *)                              fail "hook denies cargo from main checkout" "got: $dec" ;;
esac

dec=$(hook_decision "$WT_B" ".claude/bin/machine-lock run X Y -- cargo test")
case "$dec" in
    allow) pass "hook allows self-wrapping machine-lock invocation" ;;
    *)     fail "hook allows self-wrapping machine-lock invocation" "got: $dec" ;;
esac

(cd "$WT_A" && "$TEST_LOCK" release >/dev/null)

reset_locks
dec=$(hook_decision "$WT_A" "cargo test")
case "$dec" in
    deny:*machine-lock*|deny:*machine\ lock*) pass "hook denies bare cargo when no lock is held" ;;
    *)                                         fail "hook denies bare cargo when no lock is held" "got: $dec" ;;
esac


echo "==> Bug 2: atomic acquire (race-free)"

# The mkdir-based atomicity guarantees that N simultaneous acquires against
# the same lockdir cannot all succeed on the first try — exactly one mkdir
# call wins. To test this deterministically without depending on PID
# liveness races (which can steal the lock as losers retire their sleeps),
# we all share a single live owner-PID $LIVE_PID so the losers see a
# non-stale lock and return rc=1 (BLOCKED).
reset_locks
racer_pids=()
for i in 1 2 3 4 5 6 7 8; do
    (
        cd "$WT_A"
        ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire "Race-$i" "race test" >/tmp/mlfx-race-$i.out 2>&1
    ) &
    racer_pids+=($!)
done
for pid in "${racer_pids[@]}"; do wait "$pid" 2>/dev/null || true; done
winners=0
blockers=0
for i in 1 2 3 4 5 6 7 8; do
    if grep -q '^ACQUIRED' "/tmp/mlfx-race-$i.out" 2>/dev/null; then
        winners=$((winners + 1))
    elif grep -q '^BLOCKED' "/tmp/mlfx-race-$i.out" 2>/dev/null; then
        blockers=$((blockers + 1))
    fi
done
if [ "$winners" = "1" ] && [ "$blockers" = "7" ]; then
    pass "exactly one of eight racing acquires wins, others BLOCKED"
else
    fail "exactly one of eight racing acquires wins" "winners=$winners blockers=$blockers"
    for i in 1 2 3 4 5 6 7 8; do
        echo "    race-$i.out:"
        sed 's/^/      /' "/tmp/mlfx-race-$i.out"
    done
fi
(cd "$WT_A" && "$TEST_LOCK" release >/dev/null) || true
reset_locks


echo "==> Bug 3: --wait mode and serialized run"

# 3a. Non-wait acquire against a live-owner lock returns rc=1 immediately.
reset_locks
(cd "$WT_A" && ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire "Holder" "blocking" >/dev/null)

set +e
(cd "$WT_A" && "$TEST_LOCK" acquire "Hopeful" "should fail" >/dev/null 2>&1)
rc=$?
set -e
if [ "$rc" = "1" ]; then
    pass "non-wait acquire on held lock returns rc=1"
else
    fail "non-wait acquire on held lock returns rc=1" "rc=$rc"
fi

# 3b. --wait acquire succeeds after release from another process.
(
    sleep 1
    (cd "$WT_A" && "$TEST_LOCK" release >/dev/null)
) &
releaser_pid=$!

start=$(date +%s)
set +e
(cd "$WT_A" && ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire --wait --timeout=15 "Waiter" "waited test" >/dev/null 2>&1)
rc=$?
set -e
end=$(date +%s)
elapsed=$((end - start))

wait $releaser_pid 2>/dev/null || true

if [ "$rc" = "0" ] && [ "$elapsed" -ge 1 ] && [ "$elapsed" -le 9 ]; then
    pass "--wait acquires after release (elapsed=${elapsed}s)"
else
    fail "--wait acquires after release" "rc=$rc, elapsed=${elapsed}s"
fi
(cd "$WT_A" && "$TEST_LOCK" release >/dev/null) || true
reset_locks

# 3c. --wait with --timeout actually times out when nobody releases.
(cd "$WT_A" && ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire "Holder" "timeout test" >/dev/null)
start=$(date +%s)
set +e
(cd "$WT_A" && "$TEST_LOCK" acquire --wait --timeout=3 "Waiter" "timeout test" >/dev/null 2>&1)
rc=$?
set -e
end=$(date +%s)
elapsed=$((end - start))
if [ "$rc" != "0" ] && [ "$elapsed" -ge 2 ] && [ "$elapsed" -le 7 ]; then
    pass "--wait times out after --timeout seconds (elapsed=${elapsed}s)"
else
    fail "--wait times out after --timeout seconds" "rc=$rc, elapsed=${elapsed}s"
fi
(cd "$WT_A" && "$TEST_LOCK" release >/dev/null)
reset_locks

# 3d. Two concurrent `run` calls serialize: each waits for the other.
(
    cd "$WT_A"
    ("$TEST_LOCK" run "R1" "first"  -- sleep 2) >/tmp/mlfx-run1.out 2>&1 &
    p1=$!
    ("$TEST_LOCK" run "R2" "second" -- sleep 2) >/tmp/mlfx-run2.out 2>&1 &
    p2=$!
    wait "$p1" "$p2"
) || true
r1_ok=$(grep -c '^ACQUIRED' /tmp/mlfx-run1.out || true)
r2_ok=$(grep -c '^ACQUIRED' /tmp/mlfx-run2.out || true)
if [ "$r1_ok" = "1" ] && [ "$r2_ok" = "1" ]; then
    pass "concurrent run calls both eventually ACQUIRE (serialized)"
else
    fail "concurrent run calls both eventually ACQUIRE" "r1=$r1_ok r2=$r2_ok"
    echo "    run1.out: $(cat /tmp/mlfx-run1.out)"
    echo "    run2.out: $(cat /tmp/mlfx-run2.out)"
fi


echo "==> Bug 4: PID-liveness staleness"

# 4a. Lock with live PID + age > hard cap... we don't want to wait 3h.
#     Instead, verify status reports alive+not-stale for a normal live lock.
reset_locks
(cd "$WT_A" && ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire "Bencher" "simulated bench" >/dev/null)
status=$(cd "$WT_A" && "$TEST_LOCK" status 2>&1)
case "$status" in
    LOCKED*alive*) pass "live PID reports status=LOCKED with liveness=alive" ;;
    *)             fail "live PID reports status=LOCKED with liveness=alive" "status: $status" ;;
esac

# 4b. Backdate the timestamp but keep the PID alive — should still be LOCKED
#     (PID liveness beats age).
echo $(( $(date +%s) - 3600 )) > "$TEST_LOCKDIR/machine.lock.d/ts"
status=$(cd "$WT_A" && "$TEST_LOCK" status 2>&1)
case "$status" in
    LOCKED*alive*) pass "live PID + 1h age still reports LOCKED (not stale)" ;;
    *)             fail "live PID + 1h age still reports LOCKED (not stale)" "status: $status" ;;
esac

# 4c. Swap the PID to a dead one → status flips to STALE.
echo "99999999" > "$TEST_LOCKDIR/machine.lock.d/pid"
status=$(cd "$WT_A" && "$TEST_LOCK" status 2>&1)
case "$status" in
    STALE*DEAD*) pass "dead PID marks lock STALE" ;;
    *)           fail "dead PID marks lock STALE" "status: $status" ;;
esac

# 4d. A stale (dead-PID) lock is auto-cleared on next acquire.
(cd "$WT_A" && ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire "Claimer" "after dead PID" >/dev/null 2>&1)
status=$(cd "$WT_A" && "$TEST_LOCK" status 2>&1)
case "$status" in
    LOCKED*Claimer*) pass "dead-PID lock auto-clears on next acquire" ;;
    *)               fail "dead-PID lock auto-clears on next acquire" "status: $status" ;;
esac
(cd "$WT_A" && "$TEST_LOCK" release >/dev/null) || true

# 4e. Legacy flat-file lock is reported STALE + cleared on acquire.
reset_locks
printf 'OldAgent\nold-format lock\n%s\n' "$(date +%s)" > "$TEST_LOCKDIR/machine.lock"
status=$(cd "$WT_A" && "$TEST_LOCK" status 2>&1)
case "$status" in
    STALE*legacy*) pass "legacy-format lock is reported STALE" ;;
    *)             fail "legacy-format lock is reported STALE" "status: $status" ;;
esac

(cd "$WT_A" && ML_OWNER_PID=$LIVE_PID "$TEST_LOCK" acquire "New" "post-legacy" >/dev/null 2>&1)
if [ -d "$TEST_LOCKDIR/machine.lock.d" ] && [ ! -f "$TEST_LOCKDIR/machine.lock" ]; then
    pass "legacy lock is cleared on acquire"
else
    fail "legacy lock is cleared on acquire" "legacy file present: $([ -f "$TEST_LOCKDIR/machine.lock" ] && echo yes || echo no), new dir: $([ -d "$TEST_LOCKDIR/machine.lock.d" ] && echo yes || echo no)"
fi
(cd "$WT_A" && "$TEST_LOCK" release >/dev/null) || true

# 4f. Hook denies cargo from a different worktree when the lock's PID is
#     dead (treats as stale). Lock held by a process that we then kill.
reset_locks
sleep 120 &
DYING_PID=$!
(cd "$WT_A" && ML_OWNER_PID=$DYING_PID "$TEST_LOCK" acquire "Ghost" "owner will die" >/dev/null)
kill $DYING_PID 2>/dev/null || true
wait $DYING_PID 2>/dev/null || true

dec=$(hook_decision "$WT_B" "cargo test")
case "$dec" in
    deny:*[Ss]tale*) pass "hook denies cargo when lock is stale (dead PID)" ;;
    *)               fail "hook denies cargo when lock is stale (dead PID)" "got: $dec" ;;
esac


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

echo
echo "=========================================="
echo "Results: $PASS passed, $FAIL failed"
if [ "$FAIL" -gt 0 ]; then
    echo "Failed:"
    for name in "${FAIL_NAMES[@]}"; do
        echo "  - $name"
    done
    exit 1
fi
exit 0
