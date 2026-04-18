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

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    reason = "ERR2 § Phase 0: benchmarks are allowed to panic"
)]

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
    scf::{self, ScfParams, mixing::MixingMode, smearing::SmearingScheme},
    symmetry::{SpaceGroupOp, SymmetryInfo, density::symmetrize_density_g},
};
use std::collections::HashMap;

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
// WFRX Phase-1: subspace warm-start vs full dense diagonalization
//
// Measures the per-k-point eigensolve cost with a realistic "converging SCF"
// scenario: H_new is a perturbed version of H_old where the warm-start is
// H_old's eigenvectors. This matches what `diagonalize_subspace` actually
// sees in an SCF loop — the residual gate should pass once the subspace is
// close enough to invariant.
//
// Configurations span the small-molecule → production spectrum per the
// Performance Engineer playbook: n≈89 (ecut=100, sanity), n≈283 (ecut=200,
// typical), n≈893 (ecut=400, production-scale).
// ---------------------------------------------------------------------------

fn bench_wfrx_subspace(c: &mut Criterion) {
    let crystal = si_crystal();
    let pp = si_pp();
    let k = Vector3::zeros();
    let n_bands = 8;

    let mut group = c.benchmark_group("wfrx_subspace");
    group.sample_size(20);

    for &ecut in &[100.0, 200.0, 400.0] {
        let basis = BasisSet::new(&crystal.lattice, ecut);
        let n = basis.len();

        // Build a realistic H_old.
        let mut h_old = hamiltonian::build_kinetic(&basis, &k);
        let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]).unwrap();
        vnl.add_to_hamiltonian(&mut h_old, &crystal, &basis, &k);

        // Warm-start subspace: the exact lowest-n_bands eigenvectors of H_old.
        // This is the "best case" for WFRX — what the SCF sees at late
        // iterations once the density is very close to self-consistent.
        let v_prev = dense::diagonalize_lowest(&h_old, n_bands).unwrap().eigenvectors;

        // H_new: tiny Hermitian perturbation to H_old. 1e-4 on the diagonal
        // mimics a V_H / V_xc update in a converging SCF; the subspace
        // should still be close enough that the residual gate passes.
        let mut h_new = h_old.clone();
        for j in 0..n {
            h_new[(j, j)] += num_complex::Complex64::new(
                1e-4 * ((j as f64 + 1.0) * 0.123_4).sin(),
                0.0,
            );
        }

        group.bench_function(format!("full_dense_n{n}"), |b| {
            b.iter(|| {
                black_box(dense::diagonalize_lowest(black_box(&h_new), n_bands).unwrap());
            });
        });

        group.bench_function(format!("subspace_warm_n{n}"), |b| {
            b.iter(|| {
                black_box(
                    dense::diagonalize_subspace(
                        black_box(&h_new),
                        n_bands,
                        Some(black_box(&v_prev)),
                    )
                    .unwrap(),
                );
            });
        });

        // Cold-start: first-iteration path (v_prev = None). Must match
        // the full-dense baseline (no degradation) — this is the
        // correctness side of the WFRX contract.
        group.bench_function(format!("subspace_cold_n{n}"), |b| {
            b.iter(|| {
                black_box(
                    dense::diagonalize_subspace(black_box(&h_new), n_bands, None).unwrap(),
                );
            });
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

// ---------------------------------------------------------------------------
// ALOC F-5: Hamiltonian assembly — cached `faer::Mat` scratch vs per-call
// `Mat::<Complex64>::zeros(n_pw, n_pw)` allocation.
//
// The pre-ALOC-F5 driver allocated a fresh `faer::Mat<Complex64>` of size
// `n_pw × n_pw` at every SCF iteration per k-point (16·n_pw² bytes each).
// Post-ALOC-F5 the driver carries one `Mat` per k-point in
// `ScfContext::h_scratch` and fully overwrites every entry in-place via
// `fill_hamiltonian_with_v_eff`.
//
// This bench isolates the assembly cost per k-point at the three canonical
// sizes:
//   - n_pw ≈  89  (ecut = 100 eV, micro — sanity)
//   - n_pw ≈ 283  (ecut = 200 eV, medium — typical production)
//   - n_pw ≈ 893  (ecut = 400 eV, production — 16·n²=~12.7 MB/Mat)
//
// `alloc_and_fill` mimics the old path: allocate a zero Mat then run the
// legacy kinetic + `+=` assembly (equivalent to the removed
// `build_hamiltonian_with_v_eff`). `fill_into_cached` uses the new path
// with a pre-allocated buffer reused across criterion iterations. The
// difference is one allocator round-trip + one zero-write per n_pw²
// entries — the exact cost ALOC F-5 amortizes over the SCF lifetime.
// ---------------------------------------------------------------------------

fn bench_hamiltonian_assembly_aloc_f5(c: &mut Criterion) {
    use num_complex::Complex64;

    let crystal = si_crystal();
    let pp = si_pp();
    let k = Vector3::new(0.1, 0.2, 0.3); // off-Γ to exercise non-trivial kinetic
    // An FFT grid large enough to safely index every G - G' miller triple
    // at ecut = 400 eV. 32³ is the typical ecutrho = 4·ecutwfc grid for Si.
    let grid_dims = [32usize, 32, 32];

    let mut group = c.benchmark_group("aloc_f5_h_assembly");
    group.sample_size(30);

    for &ecut in &[100.0, 200.0, 400.0] {
        let basis = BasisSet::new(&crystal.lattice, ecut);
        let n = basis.len();

        // Build a realistic V_eff(G) pattern. Exact values don't matter for
        // timing — the assembly touches every entry unconditionally.
        let v_eff: Vec<Complex64> = (0..grid_dims[0] * grid_dims[1] * grid_dims[2])
            .map(|i| Complex64::new(0.001 * (i as f64 + 1.0).sin(), 0.0005 * (i as f64 + 1.0).cos()))
            .collect();

        // Pre-built VNL for post-assembly accumulation timing. The VNL cost
        // isn't ALOC F-5's target but is included so the bench is "assemble
        // the exact matrix that goes into the eigensolver," mirroring the
        // driver's call sequence.
        let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]).unwrap();

        // ---- Old path: allocate + fill + VNL each call -------------------
        group.bench_function(format!("alloc_and_fill_n{n}"), |b| {
            b.iter(|| {
                // Replicate the old `build_hamiltonian_with_v_eff` verbatim:
                // fresh zero Mat, then kinetic diagonal, then += V_eff(ΔG).
                let mut h = faer::Mat::<Complex64>::zeros(n, n);
                for (i, g) in basis.g_vectors().iter().enumerate() {
                    let ke = pwdft_rs::consts::HBAR2_OVER_2M * (k + g).norm_squared();
                    h[(i, i)] = Complex64::new(ke, 0.0);
                }
                let miller_idx = basis.miller_indices();
                for i in 0..n {
                    for j in 0..n {
                        let dn1 = i32::from(miller_idx[i][0]) - i32::from(miller_idx[j][0]);
                        let dn2 = i32::from(miller_idx[i][1]) - i32::from(miller_idx[j][1]);
                        let dn3 = i32::from(miller_idx[i][2]) - i32::from(miller_idx[j][2]);
                        // Inline miller_to_idx (it's private to scf::grid).
                        // For this FFT grid shape, wrap negative indices.
                        let wrap = |d: i32, len: usize| -> usize {
                            let l = len as i32;
                            (((d % l) + l) % l) as usize
                        };
                        let fi = wrap(dn1, grid_dims[0]);
                        let fj = wrap(dn2, grid_dims[1]);
                        let fk = wrap(dn3, grid_dims[2]);
                        let fft_idx = (fi * grid_dims[1] + fj) * grid_dims[2] + fk;
                        h[(i, j)] += v_eff[fft_idx];
                    }
                }
                vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);
                black_box(h);
            });
        });

        // ---- New path: fill into a caller-owned Mat, reused ---------------
        group.bench_function(format!("fill_into_cached_n{n}"), |b| {
            // The Mat is allocated once outside the iter loop — the
            // ALOC F-5 driver allocates once per SCF context, which
            // lives for the whole SCF run. We approximate that here with
            // a `setup + iter` pattern so the allocation isn't timed.
            let mut h = faer::Mat::<Complex64>::zeros(n, n);
            b.iter(|| {
                // In the real driver this is `fill_hamiltonian_with_v_eff`.
                // We inline the logic here so the bench is independent of
                // the crate-private helper; it is identical to the new
                // assembly path's inner loop (full overwrite, no zero-fill).
                let miller_idx = basis.miller_indices();
                let g_vectors = basis.g_vectors();
                for i in 0..n {
                    let ke_i = pwdft_rs::consts::HBAR2_OVER_2M
                        * (k + g_vectors[i]).norm_squared();
                    let mi = miller_idx[i];
                    for j in 0..n {
                        let mj = miller_idx[j];
                        let dn1 = i32::from(mi[0]) - i32::from(mj[0]);
                        let dn2 = i32::from(mi[1]) - i32::from(mj[1]);
                        let dn3 = i32::from(mi[2]) - i32::from(mj[2]);
                        let wrap = |d: i32, len: usize| -> usize {
                            let l = len as i32;
                            (((d % l) + l) % l) as usize
                        };
                        let fi = wrap(dn1, grid_dims[0]);
                        let fj = wrap(dn2, grid_dims[1]);
                        let fk = wrap(dn3, grid_dims[2]);
                        let fft_idx = (fi * grid_dims[1] + fj) * grid_dims[2] + fk;
                        let v = v_eff[fft_idx];
                        h[(i, j)] = if i == j {
                            Complex64::new(ke_i, 0.0) + v
                        } else {
                            v
                        };
                    }
                }
                vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);
                black_box(&h);
            });
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// ALOC F-5: End-to-end per-iter SCF wall-time.
//
// Runs a fixed-iteration SCF against the ALOC-F5 cache and measures the
// total time — the user-visible impact of the fix. Configurations span
// the small → medium → production spectrum required by the perf-engineer
// playbook:
//
//   - si_gamma_ecut100         : n_pw ≈  89, Γ-only           (micro)
//   - si_gamma_ecut200         : n_pw ≈ 283, Γ-only           (medium)
//   - si_2x2x2_ecut200         : n_pw ≈ 283, n_k = 4           (multi-k, medium)
//   - si_4x4x4_ecut200         : n_pw ≈ 283, n_k =10           (multi-k, production density grid)
//
// The ecut=400 configuration isn't included here — a single SCF run takes
// ~20 s/iter × 20 iters × 10 k = 4000 s, which would stall the machine-lock
// queue. The `aloc_f5_h_assembly` bench above exercises n_pw=725 on a
// per-k basis, which is where the pure allocation delta lives.
// ---------------------------------------------------------------------------

/// Reasonable SCF parameters for the ALOC F-5 end-to-end bench. Tolerance
/// is slightly loose so the `run_scf` dispatcher actually converges and
/// returns `Ok` — otherwise the bench's `.unwrap()` would fire as a
/// `ConvergenceFailure`. Once-converged, the loop exits via the normal
/// path and timing measures "full SCF to convergence," which is the
/// user-visible number.
fn build_si_scf_params(max_iter: usize) -> ScfParams {
    ScfParams {
        n_bands: 8,
        max_iter,
        conv_threshold: 1e-3,   // loose: converge in a few iters
        energy_threshold: 1e-2, // loose: energy can plateau early
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.1,
        smearing_scheme: SmearingScheme::FermiDirac,
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Plain,
        nspin: 1,
        starting_magnetization: HashMap::new(),
        ..Default::default()
    }
}

fn bench_scf_iter_end_to_end_aloc_f5(c: &mut Criterion) {
    use pwdft_rs::kpoints;

    let pp = si_pp();

    let mut group = c.benchmark_group("aloc_f5_scf_end_to_end");
    // SCF runs are long — keep criterion happy with a small sample count.
    group.sample_size(10);
    group.measurement_time(std::time::Duration::from_secs(30));

    // (label, ecut_ev, k_grid)
    let configs: &[(&str, f64, [u32; 3])] = &[
        ("si_gamma_ecut100",  100.0, [1, 1, 1]),
        ("si_gamma_ecut200",  200.0, [1, 1, 1]),
        ("si_2x2x2_ecut200",  200.0, [2, 2, 2]),
        ("si_4x4x4_ecut200",  200.0, [4, 4, 4]),
    ];

    for (label, ecut, k_grid) in configs {
        let crystal = si_crystal();
        let basis = BasisSet::new(&crystal.lattice, *ecut);
        // MPSH default: Γ-centered to match QE's `automatic 0 0 0` convention.
        let kpts = kpoints::monkhorst_pack(
            k_grid[0],
            k_grid[1],
            k_grid[2],
            kpoints::KGridShift::GammaCentered,
            &crystal.lattice,
        );
        let sym = SymmetryInfo::from_crystal(&crystal, 1e-5);
        let params = build_si_scf_params(20); // up to 20 SCF iters; converges fast

        group.bench_function(format!("{label}_5iter"), |b| {
            b.iter(|| {
                black_box(
                    scf::run_scf(
                        black_box(&crystal),
                        black_box(&basis),
                        black_box(&kpts),
                        black_box(&[&pp]),
                        black_box(&params),
                        black_box(&sym),
                    )
                    .unwrap(),
                );
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_eigensolver,
    bench_wfrx_subspace,
    bench_hamiltonian,
    bench_hamiltonian_assembly_aloc_f5,
    bench_scf_iter_end_to_end_aloc_f5,
    bench_fft,
    bench_basis,
    bench_xc_grid,
    bench_symmetry,
);
criterion_main!(benches);
