//! Criterion benchmarks for the hot computational paths.
//!
//! Run with: cargo bench
//! Results are written to target/criterion/ with HTML reports.

use criterion::{Criterion, criterion_group, criterion_main};
use nalgebra::Vector3;
use num_complex::Complex64;

use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    fft::FFT3D,
};

/// Build a standard Si FCC crystal for benchmarks.
fn si_crystal() -> Crystal {
    let a = 5.431;
    Crystal {
        lattice: Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        ),
        atoms: vec![
            Atom::new(14, [0.0, 0.0, 0.0]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    }
}

fn bench_fft_forward(c: &mut Criterion) {
    let mut group = c.benchmark_group("fft");

    for &size in &[16, 20, 24, 32] {
        let fft = FFT3D::new(size, size, size);
        let n = fft.total_size();
        let mut data: Vec<Complex64> = (0..n)
            .map(|i| Complex64::new((i as f64 * 0.1).sin(), (i as f64 * 0.2).cos()))
            .collect();

        group.bench_function(format!("forward_{size}x{size}x{size}"), |b| {
            b.iter(|| {
                fft.forward(&mut data);
            });
        });
    }

    group.finish();
}

fn bench_fft_roundtrip(c: &mut Criterion) {
    let fft = FFT3D::new(20, 20, 20);
    let n = fft.total_size();
    let mut data: Vec<Complex64> = (0..n)
        .map(|i| Complex64::new((i as f64 * 0.1).sin(), 0.0))
        .collect();

    c.bench_function("fft_roundtrip_20x20x20", |b| {
        b.iter(|| {
            fft.forward(&mut data);
            fft.inverse_normalized(&mut data);
        });
    });
}

fn bench_basis_construction(c: &mut Criterion) {
    let crystal = si_crystal();

    let mut group = c.benchmark_group("basis");
    for &ecut in &[100.0, 200.0, 400.0] {
        group.bench_function(format!("new_ecut_{ecut}"), |b| {
            b.iter(|| BasisSet::new(&crystal.lattice, ecut));
        });
    }
    group.finish();
}

fn bench_hamiltonian_build(c: &mut Criterion) {
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let n_pw = basis.len();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
    )
    .unwrap();

    let k = Vector3::new(0.0, 0.0, 0.0);

    c.bench_function(&format!("nonlocal_potential_new_{n_pw}pw"), |b| {
        b.iter(|| {
            pwdft_rs::potential::nonlocal::NonlocalPotential::new(
                &crystal,
                &basis,
                &k,
                &[&pp],
            )
        });
    });
}

criterion_group!(
    benches,
    bench_fft_forward,
    bench_fft_roundtrip,
    bench_basis_construction,
    bench_hamiltonian_build,
);
criterion_main!(benches);
