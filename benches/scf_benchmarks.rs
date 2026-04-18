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
    eigensolver::{dense, iterative},
    fft::FFT3D,
    hamiltonian,
    potential::{nonlocal::NonlocalPotential, xc},
    symmetry::{SpaceGroupOp, SymmetryInfo, density::symmetrize_density_g},
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

        // ITEV iterative partial eigensolver (lowest 8 eigenpairs only).
        // On real SCF Hamiltonians with near-degenerate eigenvalues, faer
        // 0.24's `iterate_lanczos` can spin in its inner reorthogonalization
        // loop indefinitely (see the `iterative` module doc comment for
        // details). The bench is therefore disabled on sizes that are
        // known to trigger this upstream issue until the faer fix lands.
        let _ = &iterative::DEFAULT_TOL; // keep import referenced
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

        // SCF-iteration workload: ~20 FFT calls per SCF iteration on the
        // charge-density grid (FFTB proposal). Measures buffer-reuse impact
        // on the hot SCF path; with per-call Array3 allocation this would
        // dominate the per-iteration heap traffic at large grids.
        group.bench_function(format!("scf_iter_20x_{size}x{size}x{size}"), |b| {
            b.iter(|| {
                for _ in 0..20 {
                    fft.forward(black_box(&mut data));
                    fft.inverse_normalized(black_box(&mut data));
                }
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

// ---------------------------------------------------------------------------
// LDA XC grid evaluation (per-iteration hot path)
//
// Sizes span the sequential→parallel crossover region:
//   256   — tiny molecule / sanity
//   512   — threshold candidate
//   4_096 — 16^3 grid (small SCF)
//   32_768 — 32^3 grid (typical production)
//   262_144 — 64^3 grid (large)
// ---------------------------------------------------------------------------

fn make_rho(n: usize) -> Vec<f64> {
    // Positive, physically plausible densities. Span low (near RHO_FLOOR) and
    // higher values so both the rs>=1 and rs<1 branches of PZ are exercised.
    (0..n)
        .map(|i| 0.001 + (i as f64 / n as f64) * 0.5)
        .collect()
}

fn bench_xc_grid(c: &mut Criterion) {
    let mut group = c.benchmark_group("xc_grid");

    for &n in &[256_usize, 512, 4_096, 16_384, 32_768, 262_144] {
        let rho_r = make_rho(n);

        group.bench_function(format!("lda_xc_grid_n{n}"), |b| {
            b.iter(|| black_box(xc::lda_xc_grid(black_box(&rho_r))));
        });

        let rho_up: Vec<f64> = rho_r.iter().map(|&r| 0.6 * r).collect();
        let rho_down: Vec<f64> = rho_r.iter().map(|&r| 0.4 * r).collect();

        group.bench_function(format!("lda_xc_spin_grid_n{n}"), |b| {
            b.iter(|| {
                black_box(xc::lda_xc_spin_grid(
                    black_box(&rho_up),
                    black_box(&rho_down),
                ))
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Symmetry: G-space density symmetrization
// ---------------------------------------------------------------------------

/// Build a `SymmetryInfo` with the first `k` operations from the Fd-3m group
/// (48 ops for `si_crystal()`). Used to vary `N_ops` in the symmetrizer
/// bench without changing the crystal. `k = 1` yields identity-only (which
/// the routine short-circuits), `k = 8` a sub-group subset, `k = 48` the
/// full group.
fn sym_subset(crystal: &Crystal, n_ops: usize) -> SymmetryInfo {
    let full = SymmetryInfo::from_crystal(crystal, 1e-5);
    let k = n_ops.min(full.operations.len()).max(1);
    let ops: Vec<SpaceGroupOp> = full.operations.into_iter().take(k).collect();
    SymmetryInfo {
        n_ops: ops.len(),
        operations: ops,
        has_inversion: false,
        has_time_reversal: false,
        tolerance: 1e-5,
    }
}

fn make_band_limited_rho(dims: [usize; 3]) -> Vec<f64> {
    let [nx, ny, nz] = dims;
    let n = nx * ny * nz;
    let mut rho = vec![0.0_f64; n];
    // A few low-order cosine modes so any reasonable FFT grid (n ≥ 18)
    // safely contains the rotated orbits for cubic rotations.
    let modes: [([i32; 3], f64); 5] = [
        ([0, 0, 0], 1.0),
        ([1, 0, 0], 0.17),
        ([0, 1, 0], 0.11),
        ([1, 1, 0], 0.07),
        ([1, 1, 1], 0.04),
    ];
    for ix in 0..nx {
        let fx = ix as f64 / nx as f64;
        for iy in 0..ny {
            let fy = iy as f64 / ny as f64;
            for iz in 0..nz {
                let fz = iz as f64 / nz as f64;
                let mut s = 0.0;
                for (k, amp) in &modes {
                    s += amp
                        * (std::f64::consts::TAU
                            * (k[0] as f64 * fx + k[1] as f64 * fy + k[2] as f64 * fz))
                            .cos();
                }
                rho[ix * ny * nz + iy * nz + iz] = s.abs() + 0.5;
            }
        }
    }
    rho
}

fn bench_symmetry(c: &mut Criterion) {
    let crystal = si_crystal();

    let mut group = c.benchmark_group("symmetry");
    // FFT + O(N_grid · N_ops) per call; 72³·48 dominates run-time
    // so cap samples to keep wall-time reasonable.
    group.sample_size(20);

    for &n in &[18_usize, 36, 72] {
        let dims = [n, n, n];
        let rho0 = make_band_limited_rho(dims);

        for &n_ops in &[1_usize, 8, 48] {
            let sym = sym_subset(&crystal, n_ops);
            // Reuse one FFT3D per configuration (same convention as the
            // SCF loop — FFT plans are hoisted out).
            let mut fft = FFT3D::new(dims[0], dims[1], dims[2]);
            group.bench_function(
                format!("symmetrize_density_g_n{n}_ops{n_ops}"),
                |b| {
                    b.iter_with_setup(
                        || rho0.clone(),
                        |mut rho| {
                            symmetrize_density_g(
                                black_box(&mut rho),
                                black_box(dims),
                                black_box(&mut fft),
                                black_box(&sym),
                            );
                            rho
                        },
                    );
                },
            );
        }
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_eigensolver,
    bench_hamiltonian,
    bench_fft,
    bench_basis,
    bench_xc_grid,
    bench_symmetry,
);
criterion_main!(benches);
