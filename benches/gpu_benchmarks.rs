//! GPU vs CPU benchmarks for grid-level operations.
//!
//! Run with: cargo bench --features gpu --bench gpu_benchmarks
//!
//! ## Variants
//!
//! - `cpu_<n>` — reference CPU scalar/rayon path.
//! - `gpu_fresh_<n>` — GPU kernel with per-call `create_buffer` / `create_bind_group`
//!   (the fallback when `prepare_buffers` has not been invoked). Serves as the
//!   "before" baseline for GOPT PR-B.
//! - `gpu_pooled_<n>` — GPU kernel with the `BufferPool` warm: persistent
//!   storage buffers, cached bind groups, uniform updates via `write_buffer`.
//!   This is the production path the SCF driver uses.
//!
//! Grid sizes 32³ (32_768), 64³ (262_144), 128³ (2_097_152) bracket realistic
//! SCF FFT grids — `si_scf.yaml` runs at 16³ (too small to bench GPU), and the
//! converged production runs (GGAP/QE validation class) sit in the 32³–64³
//! range.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use num_complex::Complex64;

use pwdft_rs::gpu::GpuAccelerator;
use pwdft_rs::potential::xc;

// Grid sizes spanning realistic SCF FFT meshes. 32³ = 32_768 is the smallest
// regularly-used grid; 64³ = 262_144 is the production scale for QE-validation
// runs; 128³ = 2_097_152 is the stress point above which wgpu overhead is not
// the bottleneck (F1 chain-fusion territory, deferred to GOPT PR-C).
const SIZES: &[usize] = &[32_768, 262_144, 2_097_152];

fn make_complex_data(n: usize, seed: f64) -> Vec<Complex64> {
    (0..n)
        .map(|i| {
            Complex64::new(
                (i as f64 * seed).sin() * 0.01,
                (i as f64 * seed * 1.3).cos() * 0.01,
            )
        })
        .collect()
}

fn make_g_squared(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| if i == 0 { 0.0 } else { 0.5 + i as f64 * 0.3 })
        .collect()
}

fn bench_hartree(c: &mut Criterion) {
    let Some(gpu_fresh) = GpuAccelerator::try_new() else {
        eprintln!("No GPU available, skipping GPU benchmarks");
        return;
    };

    let fourpi_e2 = 4.0 * std::f64::consts::PI * pwdft_rs::consts::E2_COULOMB;

    let mut group = c.benchmark_group("hartree");
    for &n in SIZES {
        let rho_g = make_complex_data(n, 0.1);
        let g_squared = make_g_squared(n);

        group.bench_function(format!("cpu_{n}"), |b| {
            b.iter(|| {
                let v: Vec<Complex64> = black_box(&rho_g)
                    .iter()
                    .zip(black_box(&g_squared).iter())
                    .map(|(&rho, &g2)| {
                        if g2 > 1e-20 { rho * fourpi_e2 / g2 } else { Complex64::new(0.0, 0.0) }
                    })
                    .collect();
                black_box(v);
            });
        });

        // Fresh-alloc fallback (no prepare_buffers): measures the
        // pre-GOPT-PR-B per-call buffer + bind-group creation cost.
        group.bench_function(format!("gpu_fresh_{n}"), |b| {
            b.iter(|| {
                black_box(gpu_fresh.hartree_potential(
                    black_box(&rho_g),
                    black_box(&g_squared),
                    fourpi_e2,
                ));
            });
        });

        // Pooled path: warm up a separate accelerator so the pool is sized to n.
        let mut gpu_pooled =
            GpuAccelerator::try_new().expect("second GPU init succeeded during first");
        gpu_pooled.prepare_buffers(n, &g_squared);
        group.bench_function(format!("gpu_pooled_{n}"), |b| {
            b.iter(|| {
                black_box(gpu_pooled.hartree_potential(
                    black_box(&rho_g),
                    black_box(&g_squared),
                    fourpi_e2,
                ));
            });
        });
    }
    group.finish();
}

fn bench_v_eff(c: &mut Criterion) {
    let Some(gpu_fresh) = GpuAccelerator::try_new() else { return };

    let mut group = c.benchmark_group("v_eff_assembly");
    for &n in SIZES {
        let v_local = make_complex_data(n, 0.1);
        let v_h = make_complex_data(n, 0.2);
        let v_xc = make_complex_data(n, 0.3);
        let g_squared = make_g_squared(n);

        group.bench_function(format!("cpu_{n}"), |b| {
            b.iter(|| {
                let v: Vec<Complex64> = (0..n)
                    .map(|i| black_box(&v_local)[i] + black_box(&v_h)[i] + black_box(&v_xc)[i])
                    .collect();
                black_box(v);
            });
        });

        group.bench_function(format!("gpu_fresh_{n}"), |b| {
            b.iter(|| {
                black_box(gpu_fresh.v_eff_assembly(
                    black_box(&v_local),
                    black_box(&v_h),
                    black_box(&v_xc),
                ));
            });
        });

        let mut gpu_pooled =
            GpuAccelerator::try_new().expect("second GPU init succeeded during first");
        gpu_pooled.prepare_buffers(n, &g_squared);
        group.bench_function(format!("gpu_pooled_{n}"), |b| {
            b.iter(|| {
                black_box(gpu_pooled.v_eff_assembly(
                    black_box(&v_local),
                    black_box(&v_h),
                    black_box(&v_xc),
                ));
            });
        });
    }
    group.finish();
}

fn bench_lda_xc(c: &mut Criterion) {
    let Some(gpu_fresh) = GpuAccelerator::try_new() else { return };

    let mut group = c.benchmark_group("lda_xc");
    for &n in SIZES {
        let rho_r: Vec<f64> = (0..n)
            .map(|i| 0.01 + (i as f64 / n as f64) * 0.5)
            .collect();
        let g_squared = make_g_squared(n);

        group.bench_function(format!("cpu_{n}"), |b| {
            b.iter(|| black_box(xc::lda_xc_grid(black_box(&rho_r))));
        });

        group.bench_function(format!("gpu_fresh_{n}"), |b| {
            b.iter(|| black_box(gpu_fresh.lda_xc(black_box(&rho_r))));
        });

        let mut gpu_pooled =
            GpuAccelerator::try_new().expect("second GPU init succeeded during first");
        gpu_pooled.prepare_buffers(n, &g_squared);
        group.bench_function(format!("gpu_pooled_{n}"), |b| {
            b.iter(|| black_box(gpu_pooled.lda_xc(black_box(&rho_r))));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_hartree, bench_v_eff, bench_lda_xc);
criterion_main!(benches);
