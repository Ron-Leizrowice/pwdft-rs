//! Diagnostic test: check that the non-local potential preserves
//! cubic symmetry at the Γ point for Si FCC.

use nalgebra::Vector3;
use num_complex::Complex64;

use pwdft_rs::{
    basis::BasisSet,
    consts::HBAR2_OVER_2M,
    crystal::{Atom, Crystal, Lattice},
    potential::nonlocal::NonlocalPotential,
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

#[test]
fn test_vnl_hermitian_at_gamma() {
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let k: Vector3<f64> = Vector3::zeros();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();

    let n = basis.len();
    let mut h = faer::Mat::<Complex64>::zeros(n, n);
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]).unwrap();
    vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);

    // Check Hermiticity: H(i,j) = H(j,i)*
    let mut max_err = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let diff = (h[(i, j)] - h[(j, i)].conj()).norm();
            max_err = max_err.max(diff);
        }
    }
    assert!(
        max_err < 1e-10,
        "V_NL not Hermitian at Γ: max |H(i,j) - H(j,i)*| = {max_err:.2e}"
    );
}

#[test]
fn test_vnl_diagonal_same_for_symmetry_related_g() {
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let k: Vector3<f64> = Vector3::zeros();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();

    let n = basis.len();
    let mut h_nl = faer::Mat::<Complex64>::zeros(n, n);
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]).unwrap();
    vnl.add_to_hamiltonian(&mut h_nl, &crystal, &basis, &k);

    // At Γ, G-vectors with the same |G| should have the same diagonal V_NL
    // because the angular factor P_l(cos 0) = 1 for the self-term.
    // Group G-vectors by |G|² and check diagonal elements match within each shell.
    let g_vecs = basis.g_vectors();
    let mut shells: std::collections::BTreeMap<i64, Vec<(usize, f64)>> = std::collections::BTreeMap::new();
    for (i, g) in g_vecs.iter().enumerate() {
        // Round |G|² to avoid floating-point grouping issues
        let g2_key = (g.norm_squared() * 1e6).round() as i64;
        let diag = h_nl[(i, i)].re;
        shells.entry(g2_key).or_default().push((i, diag));
    }

    for (g2_key, members) in &shells {
        if members.len() < 2 {
            continue;
        }
        let ref_val = members[0].1;
        for &(idx, val) in &members[1..] {
            let diff = (val - ref_val).abs();
            assert!(
                diff < 1e-8,
                "V_NL diagonal differs within |G|²={}: G-idx {} has {:.10}, G-idx {} has {:.10} (diff={:.2e})",
                *g2_key as f64 / 1e6,
                members[0].0, ref_val, idx, val, diff
            );
        }
    }
}

#[test]
fn test_full_hamiltonian_degeneracy_at_gamma() {
    // Build a full Hamiltonian (kinetic + V_NL only, no SCF potentials)
    // and check that symmetry-related eigenvalues are degenerate.
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let k: Vector3<f64> = Vector3::zeros();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();

    let n = basis.len();
    let mut h = faer::Mat::<Complex64>::zeros(n, n);

    // Kinetic
    for (i, g) in basis.g_vectors().iter().enumerate() {
        h[(i, i)] = Complex64::new(HBAR2_OVER_2M * g.norm_squared(), 0.0);
    }

    // Non-local only
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]).unwrap();
    vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);

    let result = pwdft_rs::eigensolver::dense::diagonalize_lowest(&h, 8).unwrap();
    eprintln!("T + V_NL eigenvalues at Γ: {:?}", result.eigenvalues);

    // Bands 2-4 should be degenerate (p-like states, O_h symmetry)
    let spread_234 = result.eigenvalues[3] - result.eigenvalues[1];
    assert!(
        spread_234 < 0.1,
        "Bands 2-4 should be degenerate at Γ: spread = {spread_234:.4} eV, values = [{:.4}, {:.4}, {:.4}]",
        result.eigenvalues[1], result.eigenvalues[2], result.eigenvalues[3]
    );
}

#[test]
fn test_local_potential_symmetry() {
    // V_local(G) should have the full symmetry of the crystal.
    // For Si FCC, G-vectors with the same |G| that are related by cubic symmetry
    // should produce V_local with the same magnitude.
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();
    let omega = crystal.lattice.volume();

    // Compute V_local at each basis G-vector directly (not via FFT grid)
    let g_vecs = basis.g_vectors();
    let mut v_local_direct: Vec<Complex64> = Vec::new();
    for g in g_vecs {
        let g_norm = g.norm();
        let mut v = Complex64::new(0.0, 0.0);
        for atom in &crystal.atoms {
            let tau = atom.cart_position(&crystal.lattice);
            let phase = -g.dot(&tau);
            let sf = Complex64::new(phase.cos(), phase.sin());
            let v_form = pp.v_local_of_g(g_norm, omega);
            v += sf * v_form;
        }
        v_local_direct.push(v);
    }

    // Group by |G|² shell and check V_local magnitude is constant within shell
    let mut shells: std::collections::BTreeMap<i64, Vec<(usize, f64)>> = std::collections::BTreeMap::new();
    for (i, g) in g_vecs.iter().enumerate() {
        let g2_key = (g.norm_squared() * 1e6).round() as i64;
        shells.entry(g2_key).or_default().push((i, v_local_direct[i].norm()));
    }

    for (g2_key, members) in &shells {
        if members.len() < 2 { continue; }
        let ref_val = members[0].1;
        for &(idx, val) in &members[1..] {
            let diff = (val - ref_val).abs();
            if diff > 1e-6 {
                let g = g_vecs[idx];
                let g0 = g_vecs[members[0].0];
                eprintln!(
                    "V_local magnitude differs in |G|²={:.4} shell: G={:.4},{:.4},{:.4} → {:.8}, G={:.4},{:.4},{:.4} → {:.8} (diff={:.2e})",
                    *g2_key as f64 / 1e6, g0.x, g0.y, g0.z, ref_val, g.x, g.y, g.z, val, diff
                );
            }
            assert!(
                diff < 1e-4,
                "|V_local| differs within shell: {ref_val:.8} vs {val:.8} (diff={diff:.2e})"
            );
        }
    }
}

#[test]
fn test_kinetic_plus_vlocal_degeneracy() {
    // Build H = T + V_local (no V_H, V_xc, or V_NL) and check degeneracy at Γ.
    // If this is broken, the issue is in V_local or the FFT grid lookup.
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let _k: Vector3<f64> = Vector3::zeros();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();
    let omega = crystal.lattice.volume();

    let n = basis.len();
    let g_vecs = basis.g_vectors();
    let _miller = basis.miller_indices();
    let mut h = faer::Mat::<Complex64>::zeros(n, n);

    // Kinetic
    for (i, g) in g_vecs.iter().enumerate() {
        h[(i, i)] = Complex64::new(HBAR2_OVER_2M * g.norm_squared(), 0.0);
    }

    // V_local via DIRECT computation (bypassing FFT grid) — gold standard
    for i in 0..n {
        for j in 0..n {
            let g_diff = g_vecs[i] - g_vecs[j];
            let g_norm = g_diff.norm();
            let mut v = Complex64::new(0.0, 0.0);
            for atom in &crystal.atoms {
                let tau = atom.cart_position(&crystal.lattice);
                let phase = -g_diff.dot(&tau);
                let sf = Complex64::new(phase.cos(), phase.sin());
                let v_form = pp.v_local_of_g(g_norm, omega);
                v += sf * v_form;
            }
            h[(i, j)] += v;
        }
    }

    let result = pwdft_rs::eigensolver::dense::diagonalize_lowest(&h, 8).unwrap();
    eprintln!("T + V_local (direct) eigenvalues at Γ: {:?}", result.eigenvalues);

    // With only T + V_local, the 3-fold degeneracy is bands 3-5 (p-like),
    // not bands 2-4 (band 2 is a separate s-like state pushed down by V_local).
    let spread_345 = result.eigenvalues[4] - result.eigenvalues[2];
    eprintln!("Bands 3-5 spread (direct V_local): {spread_345:.6} eV");
    assert!(
        spread_345 < 0.001,
        "Bands 3-5 not degenerate with direct V_local: spread={spread_345:.6}"
    );
}

#[test]
fn test_kinetic_plus_vlocal_via_fft_grid() {
    // Same as above but using the FFT grid lookup for V_local(G-G').
    // If this breaks degeneracy while direct doesn't, the FFT grid mapping is buggy.
    use pwdft_rs::fft::fft_grid_size;

    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let _k: Vector3<f64> = Vector3::zeros();
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
    )
    .unwrap();
    let omega = crystal.lattice.volume();

    let n = basis.len();
    let g_vecs = basis.g_vectors();
    let miller = basis.miller_indices();

    // Build FFT grid (same as SCF would)
    let n_max: Vec<i32> = (0..3)
        .map(|dim| miller.iter().map(|m| m[dim].abs()).max().unwrap_or(0))
        .collect();
    let scale = 2i32; // 4× ecutrho → 2× G_max
    let grid_dims = [
        fft_grid_size(scale * n_max[0]),
        fft_grid_size(scale * n_max[1]),
        fft_grid_size(scale * n_max[2]),
    ];

    let recip = crystal.lattice.reciprocal();

    // Precompute V_local on full FFT grid
    let n_grid = grid_dims[0] * grid_dims[1] * grid_dims[2];
    let mut v_local_fft = vec![Complex64::new(0.0, 0.0); n_grid];
    let [nx, ny, nz] = grid_dims;
    for (idx, v_local_val) in v_local_fft.iter_mut().enumerate() {
        let i1 = idx / (ny * nz);
        let i2 = (idx / nz) % ny;
        let i3 = idx % nz;
        let n1 = if i1 > nx / 2 { i1 as i32 - nx as i32 } else { i1 as i32 };
        let n2 = if i2 > ny / 2 { i2 as i32 - ny as i32 } else { i2 as i32 };
        let n3 = if i3 > nz / 2 { i3 as i32 - nz as i32 } else { i3 as i32 };
        let g = n1 as f64 * recip.a + n2 as f64 * recip.b + n3 as f64 * recip.c;
        let g_norm = g.norm();

        for atom in &crystal.atoms {
            let tau = atom.cart_position(&crystal.lattice);
            let phase = -g.dot(&tau);
            let sf = Complex64::new(phase.cos(), phase.sin());
            let v_form = pp.v_local_of_g(g_norm, omega);
            *v_local_val += sf * v_form;
        }
    }

    // Build H using FFT grid lookup
    let mut h = faer::Mat::<Complex64>::zeros(n, n);
    for (i, g) in g_vecs.iter().enumerate() {
        h[(i, i)] = Complex64::new(HBAR2_OVER_2M * g.norm_squared(), 0.0);
    }

    for i in 0..n {
        for j in 0..n {
            let dn1 = miller[i][0] - miller[j][0];
            let dn2 = miller[i][1] - miller[j][1];
            let dn3 = miller[i][2] - miller[j][2];
            let [nx, ny, nz] = grid_dims;
            let fi1 = ((dn1 % nx as i32) + nx as i32) as usize % nx;
            let fi2 = ((dn2 % ny as i32) + ny as i32) as usize % ny;
            let fi3 = ((dn3 % nz as i32) + nz as i32) as usize % nz;
            let fft_idx = fi1 * ny * nz + fi2 * nz + fi3;
            h[(i, j)] += v_local_fft[fft_idx];
        }
    }

    let result = pwdft_rs::eigensolver::dense::diagonalize_lowest(&h, 8).unwrap();
    eprintln!("T + V_local (FFT grid) eigenvalues at Γ: {:?}", result.eigenvalues);

    let spread_345 = result.eigenvalues[4] - result.eigenvalues[2];
    eprintln!("Bands 3-5 spread (FFT V_local): {spread_345:.6} eV");
    assert!(
        spread_345 < 0.001,
        "Bands 3-5 not degenerate with FFT V_local: spread={spread_345:.6}"
    );
}
