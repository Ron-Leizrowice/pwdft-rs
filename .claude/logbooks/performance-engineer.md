# Performance Engineer Logbook

Entries: date, measurements (actual numbers), bottleneck findings, proposals assessed. Always include hardware context.

## 2026-04-16 — Orientation (no benchmarks run yet)

**Proposal audit (code inspection only, not profiled):**
- FFTB: well-motivated, FFT allocates fresh Array3 every call. Impact unknown without profiling.
- XCPR: XC parallelization straightforward but may not help at small grid sizes (4096 points). Spin-channel join is a clean win.
- WFRX: theoretically sound but 25x claim is optimistic. Real gain likely 5-15% net.
- DVSN: strategically correct long-term target, not highest near-term impact.

**Gap found:** No proposal for allocator selection (jemalloc/mimalloc). macOS libmalloc underperforms under multithreaded pressure. Trivial integration via feature flag. Worth proposing.

**TODO next session:** Run `cargo bench`, profile SCF with Instruments.app, establish actual bottleneck distribution before proposing anything.
