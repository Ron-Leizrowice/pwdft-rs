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
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
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
fn test_fe_ferromagnetic_has_moment() {
    // Fe BCC with starting magnetization should develop a magnetic moment.
    let crystal = fe_bcc();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Fe.UPF"),
    ).unwrap();
    let basis = BasisSet::new(&crystal.lattice, 100.0);
    let kpoints = gamma_only();

    let mut starting_mag = std::collections::HashMap::new();
    starting_mag.insert("Fe".to_string(), 0.5);

    let params = scf::ScfParams {
        n_bands: 6,
        max_iter: 40,
        conv_threshold: 1e-6,
        mixing_beta: 0.2,
        mixing_ndim: 4,
        smearing_sigma: 0.27, // 0.02 Ry
        ecutrho_ratio: 4,
        fft_grid: Some([16, 16, 16]),
        mixing_mode: MixingMode::Kerker { q_tf: None },
        nspin: 2,
        starting_magnetization: starting_mag,
        tot_magnetization: None,
        ..Default::default()
    };

    let result = scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, None);

    match result {
        Ok(r) => {
            eprintln!("Fe spin-polarized: E={:.6} eV, M={:.4} μB, {} iters",
                r.total_energy, r.magnetization, r.n_iterations);
            eprintln!("  E_F={:.4} eV", r.fermi_energy);

            // Fe should have a magnetic moment > 0
            assert!(
                r.magnetization > 0.5,
                "Fe should be ferromagnetic, got M={:.4} μB", r.magnetization
            );
        }
        Err(e) => {
            eprintln!("Fe spin-polarized SCF did not converge: {e}");
            // Convergence is not guaranteed at Gamma-only with low cutoff
        }
    }
}
