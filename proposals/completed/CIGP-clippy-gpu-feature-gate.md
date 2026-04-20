---
id: CIGP
status: completed
priority: low
complexity: small
risk: low
depends_on: []
blocks: []
---

# CIGP: Document `--features gpu` in clippy CI gate

## Origin

QLN2 (PR #27, 2026-04-17) Code Reviewer flagged: the existing CLIP-era clippy gate runs `cargo clippy --all-targets` without `--features gpu`, so the GPU test binary is never linted. This is exactly how the `uninlined_format_args` violation in `tests/gpu_consistency.rs:312` slipped past CLIP and lived until QLN2 caught it.

QLN2 itself fixed the live offender. CIGP closes the process gap so future GPU-only lint hits don't silently accumulate.

## Implementation

Two paths (pick the cheaper):

1. **Documentation only.** Update `CLAUDE.md § Code Quality` to require both:

   ```bash
   cargo clippy -q --all-targets
   cargo clippy -q --all-targets --features gpu
   ```

   before declaring clippy clean. Costs nothing but relies on agents reading the doc.

2. **Hook enforcement.** Add a small wrapper script `.claude/bin/clippy-check.sh` that runs both invocations and exits non-zero on any warning. Make agent definitions reference it instead of the bare `cargo clippy` command. Stronger guarantee, slightly more infrastructure.

Recommended: option 1 first (immediate, free), reconsider hook enforcement only if it slips again.

## Verification

After option 1: a manual reread of `CLAUDE.md` confirms the dual-invocation requirement is documented in the Code Quality section.

## Notes

This proposal is process / doc, not code. Trivial — Technical Writer or Engineering Manager can land it in <15 min.
