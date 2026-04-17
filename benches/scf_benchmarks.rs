//! Criterion benchmarks for the hot computational paths.
//!
//! Run with: cargo bench --bench scf_benchmarks
//! Results are written to target/criterion/ with HTML reports.
//!
//! Benchmarks cover multiple problem sizes via energy cutoff variation:
//!   ecut=100 → n_pw≈59,  small molecule
//!   ecut=200 → n_pw≈283, typical production
//!   ecut=400 → n_pw≈893, high-accuracy
//!   ecut=600 → n_pw≈1639, stress test

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use nalgebra::Vector3;
use num_complex::Complex64;

use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    eigensolver::dense,
    fft::FFT3D,
    hamiltonian,
    potential::nonlocal::NonlocalPotential,
};

/// Si FCC crystal (2 atoms, diamond structure).
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

fn si_pp() -> pwdft_rs::pseudopotential::PseudopotentialData {
    pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap()
}

// ---------------------------------------------------------------------------
// Eigensolver: the SCF bottleneck
// ---------------------------------------------------------------------------

fn bench_eigensolver(c: &mut Criterion) {
    let crystal = si_crystal();
    let pp = si_pp();
    let k = Vector3::zeros();

    let mut group = c.benchmark_group("eigensolver");
    group.sample_size(20); // eigensolves are slow at large n

    // ecut=600 (n≈1363) crashes Accelerate in release mode due to libc++ TMO bug
    for &ecut in &[100.0, 200.0, 400.0] {
        let basis = BasisSet::new(&crystal.lattice, ecut);
        let n = basis.len();

        // Build a realistic Hamiltonian (kinetic + nonlocal, not just diagonal)
        let mut h = hamiltonian::build_kinetic(&basis, &k);
        let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]).unwrap();
        vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);

        group.bench_function(format!("faer_eigen_n{n}"), |b| {
            b.iter(|| black_box(dense::diagonalize_hermitian(black_box(&h)).unwrap()));
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Hamiltonian construction (kinetic + V_NL)
// ---------------------------------------------------------------------------

fn bench_hamiltonian(c: &mut Criterion) {
    let crystal = si_crystal();
    let pp = si_pp();
    let k_gamma = Vector3::zeros();
    let k_offgamma = Vector3::new(0.1, 0.2, 0.3);

    let mut group = c.benchmark_group("hamiltonian");

    for &ecut in &[100.0, 200.0, 400.0] {
        let basis = BasisSet::new(&crystal.lattice, ecut);
        let n = basis.len();

        group.bench_function(format!("kinetic_n{n}"), |b| {
            b.iter(|| black_box(hamiltonian::build_kinetic(&basis, &k_gamma)));
        });

        group.bench_function(format!("vnl_new_n{n}"), |b| {
            b.iter(|| {
                black_box(NonlocalPotential::new(&crystal, &basis, &k_offgamma, &[&pp]).unwrap());
            });
        });

        let h = hamiltonian::build_kinetic(&basis, &k_gamma);
        let vnl = NonlocalPotential::new(&crystal, &basis, &k_gamma, &[&pp]).unwrap();

        group.bench_function(format!("vnl_apply_n{n}"), |b| {
            b.iter(|| {
                let mut h_copy = h.clone();
                vnl.add_to_hamiltonian(&mut h_copy, &crystal, &basis, &k_gamma);
                black_box(h_copy);
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// FFT at multiple grid sizes
// ---------------------------------------------------------------------------

fn bench_fft(c: &mut Criterion) {
    let mut group = c.benchmark_group("fft");

    for &size in &[16, 20, 24, 32, 48] {
        let mut fft = FFT3D::new(size, size, size);
        let n = fft.total_size();
        let mut data: Vec<Complex64> = (0..n)
            .map(|i| Complex64::new((i as f64 * 0.1).sin(), (i as f64 * 0.2).cos()))
            .collect();

        group.bench_function(format!("forward_{size}x{size}x{size}"), |b| {
            b.iter(|| fft.forward(&mut data));
        });

        group.bench_function(format!("roundtrip_{size}x{size}x{size}"), |b| {
            b.iter(|| {
                fft.forward(&mut data);
                fft.inverse_normalized(&mut data);
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Basis construction
// ---------------------------------------------------------------------------

fn bench_basis(c: &mut Criterion) {
    let crystal = si_crystal();
    let mut group = c.benchmark_group("basis");

    for &ecut in &[100.0, 200.0, 400.0, 600.0] {
        group.bench_function(format!("new_ecut_{ecut}"), |b| {
            b.iter(|| black_box(BasisSet::new(&crystal.lattice, ecut)));
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_eigensolver,
    bench_hamiltonian,
    bench_fft,
    bench_basis,
);
criterion_main!(benches);
