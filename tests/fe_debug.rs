//! Diagnostic tests for the Fe BCC 210 eV energy discrepancy.
//! Compares each Hamiltonian component against QE individually.

use nalgebra::Vector3;
use num_complex::Complex64;
use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    eigensolver::dense,
    hamiltonian,
    potential::nonlocal::NonlocalPotential,
};

const RY_TO_EV: f64 = 13.605693122994;

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

fn fe_pp() -> pwdft_rs::pseudopotential::PseudopotentialData {
    pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Fe.UPF"),
    )
    .unwrap()
}

#[test]
fn test_fe_basis_size() {
    let crystal = fe_bcc();
    let ecut = 15.0 * RY_TO_EV; // 15 Ry
    let basis = BasisSet::new(&crystal.lattice, ecut);
    eprintln!("Fe BCC basis: {} PWs at ecut={:.1} eV ({:.1} Ry)",
        basis.len(), ecut, 15.0);
    // QE reports 79 PWs at Gamma for 15 Ry
    eprintln!("QE reports 79 PWs at Gamma");
    // Our basis should be similar (might differ by a few due to cutoff boundary)
    assert!(
        (basis.len() as i32 - 79).unsigned_abs() < 10,
        "Basis size {} differs significantly from QE's 79",
        basis.len()
    );
}

#[test]
fn test_fe_kinetic_eigenvalues() {
    // Free-electron (kinetic-only) eigenvalues must match exactly
    let crystal = fe_bcc();
    let ecut = 15.0 * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut);
    let k = Vector3::zeros();

    let h = hamiltonian::build_kinetic(&basis, &k);
    let result = dense::diagonalize_lowest(&h, 8);

    eprintln!("Fe kinetic-only eigenvalues at Gamma:");
    for (i, &e) in result.eigenvalues.iter().enumerate() {
        eprintln!("  band {i}: {e:.4} eV");
    }

    // Band 0 should be 0 (G=0)
    assert!(result.eigenvalues[0].abs() < 1e-8, "Band 0 should be ~0 eV");
}

#[test]
fn test_fe_energy_decomposition() {
    // Compare individual energy components against QE reference.
    // QE at 15 Ry:
    //   one-electron: 280.15 eV  (band energy = kinetic + V_local + V_NL from eigenvalues)
    //   hartree:        3.35 eV
    //   xc:          -299.90 eV
    //   ewald:       -584.29 eV
    //   total:       -600.93 eV (= one-electron + hartree + xc + ewald)
    //
    // QE's "total" = E_band - E_H + E_xc - E_vxc + E_ewald + V_local(G=0)*N_el
    // But QE's "one-electron" = E_band (from eigenvalues, which EXCLUDE V_local(G=0))
    // And "xc" = E_xc - E_vxc (double-counting corrected)
    //
    // So: total = one_electron + hartree + xc + ewald

    let crystal = fe_bcc();
    let pp = fe_pp();
    let _omega = crystal.lattice.volume().abs();

    // Ewald
    let e_ewald = pwdft_rs::ewald::ewald_energy(&crystal, &[&pp]);
    let qe_ewald = -42.94465594 * RY_TO_EV;
    eprintln!("Ewald:  ours={e_ewald:.4}  QE={qe_ewald:.4}  diff={:.4}", e_ewald - qe_ewald);
}

#[test]
fn test_fe_ewald_energy() {
    let crystal = fe_bcc();
    let pp = fe_pp();
    let e_ewald = pwdft_rs::ewald::ewald_energy(&crystal, &[&pp]);
    let qe_ewald = -42.94465594 * RY_TO_EV;

    eprintln!("Fe Ewald energy: {e_ewald:.6} eV");
    eprintln!("QE Ewald energy: {qe_ewald:.6} eV");
    let diff = (e_ewald - qe_ewald).abs();
    eprintln!("Diff: {diff:.6} eV");

    assert!(
        diff < 1.0,
        "Ewald energy {e_ewald:.4} eV differs from QE {qe_ewald:.4} eV by {diff:.4}"
    );
}

#[test]
fn test_fe_v_local_at_g0() {
    // V_local(G=0) is a key diagnostic — it's the average potential
    let crystal = fe_bcc();
    let pp = fe_pp();
    let omega = crystal.lattice.volume().abs();

    let v_g0 = pp.v_local_of_g(0.0, omega);
    eprintln!("Fe V_local(G=0) = {v_g0:.6} eV");
    eprintln!("Fe cell volume = {omega:.6} A^3");

    // This should be a large negative number (attractive nuclear potential)
    assert!(v_g0.is_finite(), "V_local(G=0) is not finite");
}

#[test]
fn test_fe_vnl_diagonal_at_gamma() {
    // Non-local potential matrix elements at Gamma
    let crystal = fe_bcc();
    let pp = fe_pp();
    let ecut = 15.0 * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut);
    let k = Vector3::zeros();

    let n = basis.len();
    let mut h_nl = faer::Mat::<Complex64>::zeros(n, n);
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]);
    vnl.add_to_hamiltonian(&mut h_nl, &crystal, &basis, &k);

    // Print diagonal V_NL elements (first 8)
    eprintln!("Fe V_NL diagonal at Gamma (first 8):");
    for i in 0..8.min(n) {
        eprintln!("  G={i}: V_NL = {:.6} + {:.6}i eV", h_nl[(i, i)].re, h_nl[(i, i)].im);
    }

    // V_NL should be Hermitian
    let mut max_asym = 0.0_f64;
    for i in 0..n.min(20) {
        for j in 0..n.min(20) {
            let diff = (h_nl[(i, j)] - h_nl[(j, i)].conj()).norm();
            max_asym = max_asym.max(diff);
        }
    }
    eprintln!("Max Hermiticity violation: {max_asym:.2e}");
    assert!(max_asym < 1e-10, "V_NL not Hermitian: max_asym={max_asym:.2e}");

    // Check the trace (sum of diagonal) as a sanity check
    let trace: f64 = (0..n).map(|i| h_nl[(i, i)].re).sum();
    eprintln!("V_NL trace = {trace:.6} eV");
}

#[test]
fn test_fe_full_hamiltonian_eigenvalues() {
    // Full H = kinetic + V_local + V_NL at Gamma
    // Compare against QE: 5.16  26.26  26.26  27.15  27.15  27.15  40.20  40.20
    let crystal = fe_bcc();
    let pp = fe_pp();
    let ecut = 15.0 * RY_TO_EV;
    let basis = BasisSet::new(&crystal.lattice, ecut);
    let k = Vector3::zeros();
    let omega = crystal.lattice.volume().abs();

    // Step 1: Kinetic only
    let h_kin = hamiltonian::build_kinetic(&basis, &k);
    let eig_kin = dense::diagonalize_lowest(&h_kin, 8);
    eprintln!("Kinetic-only eigenvalues: {:?}", eig_kin.eigenvalues);

    // Step 2: Kinetic + V_local (via FFT grid)
    // Build V_local on a small FFT grid
    let n = basis.len();
    let fft_dims = [16, 16, 16];
    let n_grid = fft_dims[0] * fft_dims[1] * fft_dims[2];
    let recip = crystal.lattice.reciprocal();

    let mut v_local_fft = vec![Complex64::new(0.0, 0.0); n_grid];
    for (idx, v_local_val) in v_local_fft.iter_mut().enumerate() {
        let [nx, ny, nz] = fft_dims;
        let i1 = idx / (ny * nz);
        let i2 = (idx / nz) % ny;
        let i3 = idx % nz;
        let n1 = if i1 > nx / 2 { i1 as i32 - nx as i32 } else { i1 as i32 };
        let n2 = if i2 > ny / 2 { i2 as i32 - ny as i32 } else { i2 as i32 };
        let n3 = if i3 > nz / 2 { i3 as i32 - nz as i32 } else { i3 as i32 };
        let g = n1 as f64 * recip.a + n2 as f64 * recip.b + n3 as f64 * recip.c;
        let g_norm = g.norm();

        let tau = crystal.atoms[0].cart_position(&crystal.lattice);
        let phase = -g.dot(&tau);
        let sf = Complex64::new(phase.cos(), phase.sin());
        *v_local_val = sf * pp.v_local_of_g(g_norm, omega);
    }

    // Build H = T + V_local
    let miller = basis.miller_indices();
    let mut h_loc = h_kin.clone();
    for i in 0..n {
        for j in 0..n {
            let dn = [
                miller[i][0] - miller[j][0],
                miller[i][1] - miller[j][1],
                miller[i][2] - miller[j][2],
            ];
            let [nx, ny, nz] = fft_dims;
            let i1 = ((dn[0] % nx as i32) + nx as i32) as usize % nx;
            let i2 = ((dn[1] % ny as i32) + ny as i32) as usize % ny;
            let i3 = ((dn[2] % nz as i32) + nz as i32) as usize % nz;
            let fft_idx = i1 * ny * nz + i2 * nz + i3;
            h_loc[(i, j)] += v_local_fft[fft_idx];
        }
    }
    let eig_loc = dense::diagonalize_lowest(&h_loc, 8);
    eprintln!("Kinetic + V_local eigenvalues: {:?}", eig_loc.eigenvalues);

    // Step 3: Full H = T + V_local + V_NL
    let mut h_full = h_loc.clone();
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]);
    vnl.add_to_hamiltonian(&mut h_full, &crystal, &basis, &k);
    let eig_full = dense::diagonalize_lowest(&h_full, 8);
    eprintln!("Full H eigenvalues:           {:?}", eig_full.eigenvalues);
    eprintln!("QE reference:                 [5.16, 26.26, 26.26, 27.15, 27.15, 27.15, 40.20, 40.20]");

    // Check if the pattern is qualitatively right
    let qe_eigs = [5.1621, 26.2555, 26.2555, 27.1502, 27.1502, 27.1502, 40.1971, 40.1971];
    eprintln!("\nBand-by-band comparison:");
    for (i, (&ours, &qe)) in eig_full.eigenvalues.iter().zip(qe_eigs.iter()).enumerate() {
        let diff = ours - qe;
        eprintln!("  band {i}: ours={ours:.4}  QE={qe:.4}  diff={diff:.4} eV");
    }
}
