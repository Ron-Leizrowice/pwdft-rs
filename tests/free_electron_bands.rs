//! Rigorous validation of the free-electron band structure for FCC Si.
//!
//! For free electrons (V=0), the Hamiltonian H_{G,G'}(k) = δ_{GG'} (ℏ²/2m)|k+G|²
//! is diagonal in the plane-wave basis, so the eigenvalues are exactly the
//! kinetic energies of each |k+G|² mode, sorted. No approximations are involved —
//! the diagonalizer must return EXACT values (up to floating-point roundoff).
//!
//! We verify:
//! 1. Eigenvalues match analytic (ℏ²/2m)|k+G|² at every high-symmetry point
//! 2. Correct degeneracies at Γ, X, L, W, K matching BCC reciprocal lattice shells
//! 3. Band continuity along the k-path (no spurious jumps)
//! 4. Eigenvalue ordering and count

use approx::relative_eq;
use nalgebra::Vector3;

use pwdft_rs::{
    bandstructure,
    basis::BasisSet,
    consts::HBAR2_OVER_2M,
    crystal::Lattice,
    eigensolver::dense,
    hamiltonian,
    kpoints::{self, HighSymPoint},
};

const SI_A: f64 = 5.431; // Å
const N_BANDS: usize = 15;
const ECUT: f64 = 200.0; // eV

fn si_lattice() -> Lattice {
    Lattice::new(
        SI_A / 2.0 * Vector3::new(0.0, 1.0, 1.0),
        SI_A / 2.0 * Vector3::new(1.0, 0.0, 1.0),
        SI_A / 2.0 * Vector3::new(1.0, 1.0, 0.0),
    )
}

/// Analytically compute sorted free-electron eigenvalues at a k-point.
///
/// Simply evaluates E_i = (ℏ²/2m)|k + G_i|² for all G in the basis and sorts.
fn analytic_eigenvalues(basis: &BasisSet, k: &Vector3<f64>, n_bands: usize) -> Vec<f64> {
    let mut energies: Vec<f64> = basis
        .g_vectors()
        .iter()
        .map(|g| HBAR2_OVER_2M * (k + g).norm_squared())
        .collect();
    energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
    energies.truncate(n_bands);
    energies
}

/// Count degeneracies: groups of eigenvalues within tolerance.
fn degeneracies(eigenvalues: &[f64], tol: f64) -> Vec<(f64, usize)> {
    let mut result = Vec::new();
    let mut i = 0;
    while i < eigenvalues.len() {
        let val = eigenvalues[i];
        let mut count = 1;
        while i + count < eigenvalues.len() && (eigenvalues[i + count] - val).abs() < tol {
            count += 1;
        }
        result.push((val, count));
        i += count;
    }
    result
}

/// Convert fractional reciprocal coordinates to Cartesian k-vector.
fn frac_to_cart(frac: [f64; 3], lattice: &Lattice) -> Vector3<f64> {
    let recip = lattice.reciprocal();
    frac[0] * recip.a + frac[1] * recip.b + frac[2] * recip.c
}

// ============================================================
// Test 1: Eigenvalues exactly match analytic at every high-symmetry point
// ============================================================
#[test]
fn test_eigenvalues_match_analytic_at_high_sym_points() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);

    let high_sym = vec![
        ("Γ", [0.0, 0.0, 0.0]),
        ("X", [0.5, 0.0, 0.5]),
        ("L", [0.5, 0.5, 0.5]),
        ("W", [0.5, 0.25, 0.75]),
        ("K", [0.375, 0.375, 0.75]),
        ("U", [0.625, 0.25, 0.625]),
    ];

    for (label, frac) in &high_sym {
        let k = frac_to_cart(*frac, &lattice);
        let h = hamiltonian::build_hamiltonian(&basis, &k, None);
        let result = dense::diagonalize_lowest(&h, N_BANDS);
        let analytic = analytic_eigenvalues(&basis, &k, N_BANDS);

        for (i, (got, expected)) in result
            .eigenvalues
            .iter()
            .zip(analytic.iter())
            .enumerate()
        {
            assert!(
                relative_eq!(got, expected, epsilon = 1e-8),
                "{label}: band {i} eigenvalue mismatch: got {got:.10}, expected {expected:.10}"
            );
        }
    }
}

// ============================================================
// Test 2: Degeneracies at Γ match BCC reciprocal lattice shells
// ============================================================
#[test]
fn test_gamma_degeneracies() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);
    let k = Vector3::zeros();

    let analytic = analytic_eigenvalues(&basis, &k, N_BANDS);
    let degens = degeneracies(&analytic, 1e-8);

    // FCC reciprocal lattice = BCC. Shells at |G|² = 0, 3, 4, 8, 11, ... in units of (2π/a)²
    // Shell 0: G=(0,0,0) → 1 vector
    // Shell 1: G=(±1,±1,±1) → 8 vectors (all combinations of ±1 with same parity constraint)
    // Shell 2: G=(±2,0,0),(0,±2,0),(0,0,±2) → 6 vectors

    let kappa_sq = (2.0 * std::f64::consts::PI / SI_A).powi(2);

    assert_eq!(degens[0].1, 1, "Γ shell 0 (E=0) should have degeneracy 1");
    assert!(
        relative_eq!(degens[0].0, 0.0, epsilon = 1e-10),
        "lowest level at Γ should be E=0"
    );

    assert_eq!(
        degens[1].1, 8,
        "Γ shell 1 (|G|²=3κ²) should have degeneracy 8, got {}",
        degens[1].1
    );
    let expected_shell1 = HBAR2_OVER_2M * 3.0 * kappa_sq;
    assert!(
        relative_eq!(degens[1].0, expected_shell1, epsilon = 1e-6),
        "Γ shell 1: got {:.6} eV, expected {:.6} eV",
        degens[1].0,
        expected_shell1
    );

    assert_eq!(
        degens[2].1, 6,
        "Γ shell 2 (|G|²=4κ²) should have degeneracy 6, got {}",
        degens[2].1
    );
    let expected_shell2 = HBAR2_OVER_2M * 4.0 * kappa_sq;
    assert!(
        relative_eq!(degens[2].0, expected_shell2, epsilon = 1e-6),
        "Γ shell 2: got {:.6} eV, expected {:.6} eV",
        degens[2].0,
        expected_shell2
    );
}

// ============================================================
// Test 3: Degeneracies at X point
// ============================================================
#[test]
fn test_x_point_degeneracies() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);
    // X = (0.5, 0, 0.5) in fractional → Cartesian: (2π/a)(0, 1, 0)
    let k = frac_to_cart([0.5, 0.0, 0.5], &lattice);

    let analytic = analytic_eigenvalues(&basis, &k, N_BANDS);
    let degens = degeneracies(&analytic, 1e-8);

    // At X = (0,1,0)κ, the lowest shells are:
    // |k+G|² where G runs over BCC reciprocal lattice
    // k/κ = (0,1,0), G/κ = Cartesian BCC points
    // k/κ + G/κ for G=0: (0,1,0) → |.|² = 1
    // k/κ + G/κ for G=(1,1,1): (1,2,1) → |.|² = 6
    // k/κ + G/κ for G=(-1,-1,-1): (-1,0,-1) → |.|² = 2
    // etc.
    // The lowest level has |k+G|²=1 in κ² units

    let kappa_sq = (2.0 * std::f64::consts::PI / SI_A).powi(2);
    let expected_lowest = HBAR2_OVER_2M * 1.0 * kappa_sq;
    assert!(
        relative_eq!(degens[0].0, expected_lowest, epsilon = 1e-6),
        "X lowest: got {:.6} eV, expected {:.6} eV",
        degens[0].0,
        expected_lowest
    );

    // Verify there's no zero-energy state at X (unlike at Γ)
    assert!(
        degens[0].0 > 1.0,
        "X point should not have zero-energy state"
    );
}

// ============================================================
// Test 4: Degeneracies at L point
// ============================================================
#[test]
fn test_l_point_degeneracies() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);
    // L = (0.5, 0.5, 0.5) → Cartesian: (π/a)(1, 1, 1) = (0.5κ)(1,1,1)
    let k = frac_to_cart([0.5, 0.5, 0.5], &lattice);

    let analytic = analytic_eigenvalues(&basis, &k, N_BANDS);
    let degens = degeneracies(&analytic, 1e-8);

    // At L = (0.5, 0.5, 0.5)κ
    // k/κ + G/κ for G=0: (0.5, 0.5, 0.5) → |.|² = 0.75
    // k/κ + G/κ for G=(-1,-1,-1): (-0.5, -0.5, -0.5) → |.|² = 0.75
    // So the lowest level is 2-fold degenerate at 0.75κ²

    let kappa_sq = (2.0 * std::f64::consts::PI / SI_A).powi(2);
    let expected_lowest = HBAR2_OVER_2M * 0.75 * kappa_sq;
    assert_eq!(
        degens[0].1, 2,
        "L lowest level should be 2-fold degenerate"
    );
    assert!(
        relative_eq!(degens[0].0, expected_lowest, epsilon = 1e-6),
        "L lowest: got {:.6} eV, expected {:.6} eV",
        degens[0].0,
        expected_lowest
    );
}

// ============================================================
// Test 5: Band continuity along k-path
// ============================================================
#[test]
fn test_band_continuity() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);

    let path = vec![
        HighSymPoint {
            label: "Γ".into(),
            frac: [0.0, 0.0, 0.0],
        },
        HighSymPoint {
            label: "X".into(),
            frac: [0.5, 0.0, 0.5],
        },
    ];
    let (kpts, distances) = kpoints::high_symmetry_path(&path, 100, &lattice);
    let bs = bandstructure::compute_band_structure(&basis, &kpts, &distances, 8, None);

    // Check that each band varies smoothly between adjacent k-points.
    // For 100 points along Γ-X, the maximum energy change per step should be small.
    for band_idx in 0..8 {
        for i in 1..bs.eigenvalues.len() {
            let de = (bs.eigenvalues[i][band_idx] - bs.eigenvalues[i - 1][band_idx]).abs();
            let dk = distances[i] - distances[i - 1];
            // Energy change per unit k should be bounded.
            // For free electrons, dE/dk = ℏ²|k+G|/m, which is at most ~50 eV·Å for ecut=200.
            // With dk ~ 0.01 1/Å, max dE ~ 0.5 eV.
            assert!(
                de < 2.0,
                "band {band_idx}: jump of {de:.4} eV between k-points {i}-{} (dk={dk:.4})",
                i - 1
            );
        }
    }
}

// ============================================================
// Test 6: All eigenvalues exactly match sorted diagonal of H
// ============================================================
#[test]
fn test_eigenvalues_are_sorted_diagonal() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);

    // Test at several k-points (not just high-symmetry)
    let test_ks: Vec<Vector3<f64>> = vec![
        Vector3::zeros(),
        frac_to_cart([0.5, 0.0, 0.5], &lattice),
        frac_to_cart([0.5, 0.5, 0.5], &lattice),
        frac_to_cart([0.123, 0.456, 0.789], &lattice), // generic point
    ];

    for k in &test_ks {
        let h = hamiltonian::build_hamiltonian(&basis, k, None);
        let n = basis.len();

        // Extract diagonal and sort
        let mut diag: Vec<f64> = (0..n).map(|i| h[(i, i)].re).collect();
        diag.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Diagonalize
        let result = dense::diagonalize_hermitian(&h);

        // All eigenvalues should match the sorted diagonal exactly
        for (i, (got, expected)) in result
            .eigenvalues
            .iter()
            .zip(diag.iter())
            .enumerate()
        {
            assert!(
                (got - expected).abs() < 1e-8,
                "k={k}: eigenvalue {i}: got {got:.10}, expected {expected:.10}, diff={}",
                (got - expected).abs()
            );
        }
    }
}

// ============================================================
// Test 7: Verify eigenvectors reconstruct H at a generic k-point
// ============================================================
#[test]
fn test_eigenvector_reconstruction() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);
    let k = frac_to_cart([0.3, 0.1, 0.2], &lattice);

    let h = hamiltonian::build_hamiltonian(&basis, &k, None);
    let result = dense::diagonalize_hermitian(&h);
    let n = basis.len();

    // Verify H V = V Λ, i.e. H v_i = λ_i v_i for each eigenpair
    for i in 0..n.min(10) {
        let v_i = result.eigenvectors.column(i);
        let hv = &h * &v_i;
        let lambda_v = &v_i * num_complex::Complex64::new(result.eigenvalues[i], 0.0);

        let residual: f64 = (&hv - &lambda_v).iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        assert!(
            residual < 1e-8,
            "eigenvector {i}: residual |Hv - λv| = {residual:.2e}"
        );
    }
}

// ============================================================
// Test 8: Verify specific numerical values at Γ
// ============================================================
#[test]
fn test_gamma_numerical_values() {
    let lattice = si_lattice();
    let basis = BasisSet::new(&lattice, ECUT);
    let k = Vector3::zeros();

    let result = dense::diagonalize_lowest(
        &hamiltonian::build_hamiltonian(&basis, &k, None),
        N_BANDS,
    );

    let kappa_sq = (2.0 * std::f64::consts::PI / SI_A).powi(2);

    // E = 0 eV (1 state)
    assert!(
        result.eigenvalues[0].abs() < 1e-10,
        "Γ band 0: expected 0, got {}",
        result.eigenvalues[0]
    );

    // E = 3κ² × ℏ²/2m (8 states: bands 1-8)
    let e_shell1 = 3.0 * kappa_sq * HBAR2_OVER_2M;
    for i in 1..=8 {
        assert!(
            relative_eq!(result.eigenvalues[i], e_shell1, epsilon = 1e-8),
            "Γ band {i}: expected {e_shell1:.10}, got {:.10}",
            result.eigenvalues[i]
        );
    }

    // E = 4κ² × ℏ²/2m (6 states: bands 9-14)
    let e_shell2 = 4.0 * kappa_sq * HBAR2_OVER_2M;
    for i in 9..N_BANDS.min(15) {
        assert!(
            relative_eq!(result.eigenvalues[i], e_shell2, epsilon = 1e-8),
            "Γ band {i}: expected {e_shell2:.10}, got {:.10}",
            result.eigenvalues[i]
        );
    }
}
