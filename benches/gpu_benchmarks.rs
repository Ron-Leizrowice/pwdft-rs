//! GPU vs CPU benchmarks for grid-level operations.
//!
//! Run with: cargo bench --features gpu --bench gpu_benchmarks

use criterion::{Criterion, criterion_group, criterion_main, black_box};
use num_complex::Complex64;

use pwdft_rs::gpu::GpuAccelerator;
use pwdft_rs::potential::xc;

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
    let Some(gpu) = GpuAccelerator::try_new() else {
        eprintln!("No GPU available, skipping GPU benchmarks");
        return;
    };

    let fourpi_e2 = 4.0 * std::f64::consts::PI * pwdft_rs::consts::E2_COULOMB;

    let mut group = c.benchmark_group("hartree");
    for &n in &[8_000, 64_000, 512_000] {
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

        group.bench_function(format!("gpu_{n}"), |b| {
            b.iter(|| {
                black_box(gpu.hartree_potential(black_box(&rho_g), black_box(&g_squared), fourpi_e2));
            });
        });
    }
    group.finish();
}

fn bench_v_eff(c: &mut Criterion) {
    let Some(gpu) = GpuAccelerator::try_new() else { return };

    let mut group = c.benchmark_group("v_eff_assembly");
    for &n in &[8_000, 64_000, 512_000] {
        let v_local = make_complex_data(n, 0.1);
        let v_h = make_complex_data(n, 0.2);
        let v_xc = make_complex_data(n, 0.3);

        group.bench_function(format!("cpu_{n}"), |b| {
            b.iter(|| {
                let v: Vec<Complex64> = (0..n)
                    .map(|i| black_box(&v_local)[i] + black_box(&v_h)[i] + black_box(&v_xc)[i])
                    .collect();
                black_box(v);
            });
        });

        group.bench_function(format!("gpu_{n}"), |b| {
            b.iter(|| {
                black_box(gpu.v_eff_assembly(black_box(&v_local), black_box(&v_h), black_box(&v_xc)));
            });
        });
    }
    group.finish();
}

fn bench_lda_xc(c: &mut Criterion) {
    let Some(gpu) = GpuAccelerator::try_new() else { return };

    let mut group = c.benchmark_group("lda_xc");
    for &n in &[8_000, 64_000, 512_000] {
        let rho_r: Vec<f64> = (0..n)
            .map(|i| 0.01 + (i as f64 / n as f64) * 0.5)
            .collect();

        group.bench_function(format!("cpu_{n}"), |b| {
            b.iter(|| black_box(xc::lda_xc_grid(black_box(&rho_r))));
        });

        group.bench_function(format!("gpu_{n}"), |b| {
            b.iter(|| black_box(gpu.lda_xc(black_box(&rho_r))));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_hartree, bench_v_eff, bench_lda_xc);
criterion_main!(benches);
