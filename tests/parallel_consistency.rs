//! Verify that parallel (rayon) and single-threaded execution produce
//! identical numerical results for FFT, density, and the full SCF pipeline.

use nalgebra::Vector3;
use num_complex::Complex64;

use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    fft::FFT3D,
};

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

fn make_test_data(n: usize) -> Vec<Complex64> {
    (0..n)
        .map(|i| {
            Complex64::new(
                (i as f64 * 0.123).sin(),
                (i as f64 * 0.456).cos(),
            )
        })
        .collect()
}

#[test]
fn test_fft_serial_vs_parallel() {
    let original = make_test_data(20 * 20 * 20);

    let mut data_serial = original.clone();
    let result_serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| {
            let mut fft = FFT3D::new(20, 20, 20);
            fft.forward(&mut data_serial);
            data_serial.clone()
        });

    let mut fft = FFT3D::new(20, 20, 20);
    let mut data_parallel = original;
    fft.forward(&mut data_parallel);

    for (i, (s, p)) in result_serial.iter().zip(data_parallel.iter()).enumerate() {
        assert!(
            (s - p).norm() < 1e-12,
            "FFT mismatch at index {i}: serial={s}, parallel={p}"
        );
    }
}

#[test]
fn test_fft_inverse_serial_vs_parallel() {
    let original = make_test_data(20 * 20 * 20);

    let mut data_serial = original.clone();
    let result_serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| {
            let mut fft = FFT3D::new(20, 20, 20);
            fft.inverse_normalized(&mut data_serial);
            data_serial.clone()
        });

    let mut fft = FFT3D::new(20, 20, 20);
    let mut data_parallel = original;
    fft.inverse_normalized(&mut data_parallel);

    for (i, (s, p)) in result_serial.iter().zip(data_parallel.iter()).enumerate() {
        assert!(
            (s - p).norm() < 1e-12,
            "Inverse FFT mismatch at index {i}: serial={s}, parallel={p}"
        );
    }
}

#[test]
fn test_fft_roundtrip_preserves_data() {
    let mut fft = FFT3D::new(20, 20, 20);
    let original = make_test_data(fft.total_size());

    let mut data = original.clone();
    fft.forward(&mut data);
    fft.inverse_normalized(&mut data);

    for (i, (got, want)) in data.iter().zip(original.iter()).enumerate() {
        assert!(
            (got - want).norm() < 1e-10,
            "Roundtrip failed at {i}: got={got}, want={want}"
        );
    }
}

#[test]
fn test_scf_serial_vs_parallel() {
    // Compare single-threaded vs multi-threaded SCF over a short run.
    // Uses Γ-only (1 k-point) and 5 iterations to keep debug-mode runtime
    // under 10 seconds while still exercising the full SCF pipeline.
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 100.0); // smaller basis for speed
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();

    // Γ-only: 1 k-point, fast but still exercises all code paths
    let kpoints = vec![pwdft_rs::kpoints::KPoint {
        k: Vector3::zeros(),
        weight: 1.0,
        label: None,
    }];

    let params = pwdft_rs::scf::ScfParams {
        n_bands: 4,
        max_iter: 5,
        conv_threshold: 1e-20, // Won't converge in 5 iters — that's fine
        mixing_beta: 0.3,
        mixing_ndim: 4,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        ..Default::default()
    };

    // Single-threaded: run 5 SCF iterations
    let result_serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| {
            pwdft_rs::scf::run_scf(
                &crystal, &basis, &kpoints, &[&pp], &params,
                &pwdft_rs::symmetry::SymmetryInfo::identity_only(),
            )
        });

    // Multi-threaded: same 5 iterations
    let result_parallel =
        pwdft_rs::scf::run_scf(
                &crystal, &basis, &kpoints, &[&pp], &params,
                &pwdft_rs::symmetry::SymmetryInfo::identity_only(),
            );

    // TAUD finding 5.3: pattern-match the specific ConvergenceFailure variant
    // rather than any Err. If the SCF started returning e.g. Eigensolver or
    // Gpu errors for unrelated reasons, `is_err()` would still pass and hide
    // the regression. `ConvergenceFailure` is specifically the expected
    // outcome from max_iter=5 with conv_threshold=1e-20 on Si.
    match &result_serial {
        Err(pwdft_rs::error::PwdftError::ConvergenceFailure { .. }) => {}
        Ok(_) => panic!(
            "serial 5-iter SCF unexpectedly converged — conv_threshold=1e-20 \
             should force ConvergenceFailure in 5 iters"
        ),
        Err(other) => panic!(
            "serial 5-iter SCF returned unexpected error variant: {other}. \
             Expected ConvergenceFailure."
        ),
    }
    match &result_parallel {
        Err(pwdft_rs::error::PwdftError::ConvergenceFailure { .. }) => {}
        Ok(_) => panic!(
            "parallel 5-iter SCF unexpectedly converged — conv_threshold=1e-20 \
             should force ConvergenceFailure in 5 iters"
        ),
        Err(other) => panic!(
            "parallel 5-iter SCF returned unexpected error variant: {other}. \
             Expected ConvergenceFailure."
        ),
    }

    // Now run to convergence with a small but real problem
    let params_conv = pwdft_rs::scf::ScfParams {
        n_bands: 4,
        max_iter: 40,
        conv_threshold: 1e-6,
        ..params
    };

    let result_s = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| {
            pwdft_rs::scf::run_scf(
                &crystal, &basis, &kpoints, &[&pp], &params_conv,
                &pwdft_rs::symmetry::SymmetryInfo::identity_only(),
            )
        });

    let result_p =
        pwdft_rs::scf::run_scf(
            &crystal, &basis, &kpoints, &[&pp], &params_conv,
            &pwdft_rs::symmetry::SymmetryInfo::identity_only(),
        );

    let s = result_s.expect("serial SCF must converge");
    let p = result_p.expect("parallel SCF must converge");

    // Convergence guard (TAUD finding 5.1): hitting max_iter is silent failure.
    assert!(
        s.n_iterations < params_conv.max_iter,
        "serial SCF hit max_iter={} without converging",
        params_conv.max_iter
    );
    assert!(
        p.n_iterations < params_conv.max_iter,
        "parallel SCF hit max_iter={} without converging",
        params_conv.max_iter
    );

    assert_eq!(s.eigenvalues.len(), p.eigenvalues.len());
    for (ik, (evs_s, evs_p)) in s.eigenvalues.iter().zip(p.eigenvalues.iter()).enumerate() {
        for (ib, (&es, &ep)) in evs_s.iter().zip(evs_p.iter()).enumerate() {
            assert!(
                (es - ep).abs() < 1e-6,
                "Eigenvalue mismatch at k={ik} band={ib}: serial={es:.6}, parallel={ep:.6}"
            );
        }
    }
    assert!(
        (s.total_energy - p.total_energy).abs() < 1e-4,
        "Total energy mismatch: serial={:.6}, parallel={:.6}",
        s.total_energy,
        p.total_energy
    );
}

#[test]
fn test_scf_kerker_serial_vs_parallel() {
    // Same as above but with Kerker preconditioning enabled.
    // Exercises the FFT→filter→IFFT path inside the mixer under parallelism.
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();

    let kpoints = vec![pwdft_rs::kpoints::KPoint {
        k: Vector3::zeros(),
        weight: 1.0,
        label: None,
    }];

    let params = pwdft_rs::scf::ScfParams {
        n_bands: 4,
        max_iter: 40,
        conv_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 4,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        mixing_mode: pwdft_rs::scf::mixing::MixingMode::Kerker { q_tf: None },
        ..Default::default()
    };

    let result_s = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| {
            pwdft_rs::scf::run_scf(
                &crystal, &basis, &kpoints, &[&pp], &params,
                &pwdft_rs::symmetry::SymmetryInfo::identity_only(),
            )
        });

    let result_p =
        pwdft_rs::scf::run_scf(
                &crystal, &basis, &kpoints, &[&pp], &params,
                &pwdft_rs::symmetry::SymmetryInfo::identity_only(),
            );

    let s = result_s.expect("serial Kerker SCF must converge");
    let p = result_p.expect("parallel Kerker SCF must converge");

    // Convergence guard (TAUD finding 5.1).
    assert!(
        s.n_iterations < params.max_iter,
        "serial Kerker SCF hit max_iter={} without converging",
        params.max_iter
    );
    assert!(
        p.n_iterations < params.max_iter,
        "parallel Kerker SCF hit max_iter={} without converging",
        params.max_iter
    );

    for (ik, (evs_s, evs_p)) in s.eigenvalues.iter().zip(p.eigenvalues.iter()).enumerate() {
        for (ib, (&es, &ep)) in evs_s.iter().zip(evs_p.iter()).enumerate() {
            assert!(
                (es - ep).abs() < 1e-6,
                "Kerker eigenvalue mismatch at k={ik} band={ib}: serial={es:.6}, parallel={ep:.6}"
            );
        }
    }
    assert!(
        (s.total_energy - p.total_energy).abs() < 1e-4,
        "Kerker energy mismatch: serial={:.6}, parallel={:.6}",
        s.total_energy,
        p.total_energy
    );
}
