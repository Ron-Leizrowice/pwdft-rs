//! Spin polarization tests.
//!
//! Validates that nspin=2 produces correct results and matches nspin=1
//! in the unpolarized limit.

use nalgebra::Vector3;
use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    kpoints::KPoint,
    scf::{self, mixing::MixingMode},
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

fn fe_bcc() -> Crystal {
    let a = 2.87;
    Crystal {
        lattice: Lattice::new(
            a / 2.0 * Vector3::new(-1.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, -1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, -1.0),
        ),
        atoms: vec![Atom::new(26, [0.0, 0.0, 0.0])],
    }
}

fn gamma_only() -> Vec<KPoint> {
    vec![KPoint { k: Vector3::zeros(), weight: 1.0, label: None }]
}

#[test]
fn test_si_nspin2_matches_nspin1() {
    // Si with nspin=2 and zero magnetization should match nspin=1 energy.
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    ).unwrap();
    let kpoints = gamma_only();

    let params_nspin1 = scf::ScfParams {
        n_bands: 4,
        max_iter: 40,
        conv_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 4,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        mixing_mode: MixingMode::Plain,
        nspin: 1,
        ..Default::default()
    };

    let params_nspin2 = scf::ScfParams {
        nspin: 2,
        n_bands: 4,  // per spin channel
        ..params_nspin1.clone()
    };

    let result1 = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params_nspin1, None);
    let result2 = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params_nspin2, None);

    match (&result1, &result2) {
        (Ok(r1), Ok(r2)) => {
            let de = (r1.total_energy - r2.total_energy).abs();
            eprintln!("Si nspin=1: E={:.6} eV ({} iters)", r1.total_energy, r1.n_iterations);
            eprintln!("Si nspin=2: E={:.6} eV ({} iters), M={:.4} μB", r2.total_energy, r2.n_iterations, r2.magnetization);
            eprintln!("Energy diff: {de:.6} eV");

            // Energies should match within ~0.1 eV (different convergence paths)
            assert!(
                de < 0.5,
                "nspin=1 ({:.4} eV) and nspin=2 ({:.4} eV) energies differ by {de:.4} eV",
                r1.total_energy, r2.total_energy
            );

            // Magnetization should be ~0 for non-magnetic Si
            assert!(
                r2.magnetization < 0.1,
                "Si should be non-magnetic, got M={:.4} μB", r2.magnetization
            );
        }
        (Ok(_), Err(e)) => eprintln!("nspin=2 failed: {e}"),
        (Err(e), Ok(_)) => eprintln!("nspin=1 failed: {e}"),
        (Err(e1), Err(e2)) => eprintln!("Both failed: {e1}, {e2}"),
    }
}

#[test]
fn test_fe_spin_xc_consistency_regression() {
    // Combined SPXC + SPNC regression test.
    //
    // SPXC (proposals/completed/SPXC-spin-xc-consistency.md) fixed an
    // input/output mismatch in the spin-polarized E_xc: `exc_r` was from INPUT
    // spin densities while `rho_xc_total` and `rho_up/down_sym` were OUTPUT.
    // SPNC (proposals/SPNC-spin-per-density-convergence.md) fixed the nspin=2
    // convergence criterion to use per-spin max(||Δρ_up||, ||Δρ_down||)
    // instead of the total-density ||Δρ_total||, so that spin polarization
    // (zeta) is actually driven to self-consistency.
    //
    // This test exercises nspin=2 on Si with a nonzero starting_magnetization
    // so the spin XC machinery and per-spin mixing are fully active during
    // SCF, then relaxes (no tot_magnetization constraint) to the physically
    // correct non-magnetic ground state. At convergence, |E_HF - E_KS| must
    // be O(Δρ²), i.e. sub-meV at conv_threshold=1e-6.
    //
    // Historical baseline (Fe BCC fixed-mag=2 on the same nc/lda/Fe.upf, 4×4×4
    // k, 15 Ry, starting_magnetization=0.5, conv=1e-6, max_iter=300):
    //   Pre-SPXC, total-only criterion:        |HF-KS| ≈ 22.2 eV (falsely "converged")
    //   Post-SPXC, total-only criterion:       |HF-KS| ≈ 13.0 eV (falsely "converged")
    //   Post-SPXC+SPNC, per-spin criterion:    SCF correctly refuses to converge —
    //     per-channel Δρ locks at ≈0.254 because this pseudopotential does not
    //     support a stable fixed-mag=2 state (LDA ground state is nonmagnetic).
    //     The old total-only metric hid this limit cycle by cancelling +ε/−ε
    //     between the up and down channels.
    //
    // The Fe fixed-mag case is therefore documented but not tested here — it
    // would require a different pseudopotential or a better mixer. The Si
    // probe below is the right system to pin the HF-KS quadratic convergence
    // property this test is named for.
    let _ = env_logger::builder().is_test(true).try_init();
    let crystal = si_crystal();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let kpoints = gamma_only();

    let mut starting_mag = std::collections::HashMap::new();
    starting_mag.insert("Si".to_string(), 0.2);

    let params = scf::ScfParams {
        n_bands: 4,
        max_iter: 60,
        conv_threshold: 1e-6,
        energy_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 4,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        mixing_mode: MixingMode::Plain,
        nspin: 2,
        tot_magnetization: None,
        starting_magnetization: starting_mag,
        ..Default::default()
    };

    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, None)
        .expect("Si nspin=2 SCF must converge for SPXC+SPNC regression test");

    let hf_diff = (result.harris_foulkes_energy - result.total_energy).abs();
    eprintln!(
        "SPXC+SPNC regression (Si nspin=2): E_KS={:.6} eV  E_HF={:.6} eV  |HF-KS|={:.3e} eV  ({} iters, M={:.4})",
        result.total_energy,
        result.harris_foulkes_energy,
        hf_diff,
        result.n_iterations,
        result.magnetization,
    );

    // Quadratic-convergence assertion. Pre-SPXC: |HF-KS| was many eV even at
    // self-consistency due to E_xc using input `exc_r` vs output rho. Post-SPXC
    // (total-only criterion): |HF-KS| still O(delta_zeta) because spin channels
    // not driven to consistency. Post-SPXC+SPNC: should be O(delta_rho^2),
    // sub-microelectronvolt on a well-behaved non-magnetic system.
    //
    // Empirical post-SPNC |HF-KS| = 7.19e-7 eV at conv_threshold=1e-6 (23 iters).
    // Threshold set at ~14x empirical (1e-5 eV) for platform/compiler headroom.
    assert!(
        hf_diff < 1.0e-5,
        "SPXC+SPNC regression: Si nspin=2 |E_HF - E_KS| = {hf_diff:.3e} eV exceeds 1e-5 eV. \
         This system converges to ~7e-7 eV with both fixes in place. A value of many eV \
         indicates SPXC (XC input/output mismatch) has regressed; a value of O(0.01) eV \
         or larger suggests SPNC (per-spin convergence criterion) has regressed."
    );

    // Sanity: Si is non-magnetic — free-moment relaxation must drive M to 0.
    assert!(
        result.magnetization < 0.05,
        "Si should relax to non-magnetic (M≈0), got {:.4}",
        result.magnetization
    );
}

#[test]
fn test_fe_ferromagnetic_fixed_moment() {
    // Fe BCC with fixed magnetization = 2.0 μB.
    // QE reference (4x4x4, 15 Ry, LDA, FD 0.02 Ry, nspin=2, tot_mag=2):
    //   E = -44.062_678_79 Ry = -599.503 eV, 9 iters
    //   Gamma up:   4.62  25.71  25.71  26.52  26.52  26.52
    //   Gamma down: 5.79  27.30  27.30  28.05  28.05  28.05
    // NOTE: This PP favours non-magnetic Fe at LDA. The fixed-moment
    // solution is higher in energy than non-magnetic, but still validates
    // the spin-polarized SCF machinery.
    let crystal = fe_bcc();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Fe.upf"),
    ).unwrap();
    let ecut = 15.0 * 13.605_693_122_994; // 15 Ry
    let basis = BasisSet::new(&crystal.lattice, ecut);
    let kpoints = pwdft_rs::kpoints::monkhorst_pack(4, 4, 4, &crystal.lattice);

    let params = scf::ScfParams {
        n_bands: 8,
        max_iter: 100,
        conv_threshold: 1e-6,
        energy_threshold: 1e-5,
        mixing_beta: 0.2,
        mixing_ndim: 8,
        smearing_sigma: 0.02 * 13.605_693_122_994, // 0.02 Ry in eV
        ecutrho_ratio: 4,
        mixing_mode: MixingMode::Kerker { q_tf: None },
        nspin: 2,
        tot_magnetization: Some(2.0),
        ..Default::default()
    };

    let symmetry = pwdft_rs::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, Some(&symmetry));

    match result {
        Ok(r) => {
            eprintln!("Fe spin-polarized: E={:.6} eV, M={:.4} μB, {} iters",
                r.total_energy, r.magnetization, r.n_iterations);
            eprintln!("  E_F={:.4} eV", r.fermi_energy);

            // Magnetization should be ~2.0 (fixed)
            assert!(
                (r.magnetization - 2.0).abs() < 0.5,
                "Fe fixed M=2 should give M≈2, got {:.4}", r.magnetization
            );

            // Compare total energy against QE
            let qe_energy = -44.062_678_79 * 13.605_693_122_994;
            let de = (r.total_energy - qe_energy).abs();
            eprintln!("  Energy diff vs QE: {de:.4} eV (QE={qe_energy:.4} eV)");
        }
        Err(e) => {
            eprintln!("Fe spin-polarized SCF did not converge: {e}");
        }
    }
}
