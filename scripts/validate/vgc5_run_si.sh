#!/usr/bin/env bash
# VGC5 orchestrator: compute per-component energy decomposition for Si
# and Fe and compare against QE reference.
#
# 1. Parse QE standard outputs into CSVs (vgc5_qe_si_components.csv,
#    vgc5_qe_fe_components.csv).
# 2. Run the pwdft-rs integration test that prints per-component values
#    alongside QE and pins the Rust numbers.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)

cd "$REPO"

echo "==> Parsing QE outputs"
"$HERE/vgc5_per_component.py"

echo
echo "==> Running pwdft-rs per-component test (Si + Fe)"
"$REPO/.claude/bin/machine-lock" acquire "VGC5" "cargo test --release vgc5_" >/dev/null
trap '"$REPO/.claude/bin/machine-lock" release >/dev/null 2>&1 || true' EXIT
cargo test --release --test vgc5_per_component_si -- --nocapture
