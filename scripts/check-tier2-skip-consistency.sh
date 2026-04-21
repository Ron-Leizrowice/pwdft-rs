#!/usr/bin/env bash
# Check that every SKIP-TIER2 marker in #[ignore] strings corresponds to a
# real fn declaration in the next non-attribute line. Print the enumerated
# skip list. Exit 0 on success, 1 if any marker is orphaned.
#
# Handles multi-line #[ignore = "..."] strings where the value spans
# multiple lines (backslash continuation ending with "] on the last line).
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
echo "=== SKIP-TIER2 enumerated list ==="
awk '
  /SKIP-TIER2/ { found=1; in_str=(/"\]$/ ? 0 : 1); next }
  found && in_str { in_str=(/"\]$/ ? 0 : 1); next }
  found && /^[[:space:]]*(#\[|pub[[:space:]]+)?fn[[:space:]]+[a-zA-Z_]+/ {
    match($0, /fn [a-zA-Z_]+/)
    print "  --skip " substr($0, RSTART+3, RLENGTH-3)
    found=0
  }
  found && /^[[:space:]]*(#\[|\/\/)/ { next }
  found {
    print "ERROR: SKIP-TIER2 marker not followed by fn (orphan) at: " $0
    exit 1
  }
' pwdft/pwdft-core/tests/*.rs | sort -u
echo "=== done ==="
