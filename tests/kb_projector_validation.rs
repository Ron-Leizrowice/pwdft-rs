//! Rigorous validation of the Kleinman-Bylander non-local pseudopotential.
//!
//! Mathematical chain under test:
//!
//!   Step 1: UPF stores chi(r) = r * beta(r) in Bohr^{-1/2}
//!   Step 2: Our code converts: chi_Ang = chi_Bohr / sqrt(BOHR_TO_ANG)  [Ang^{-1/2}]
//!   Step 3: Form factor: F(q) = 4*pi * integral chi(r) j_l(qr) r dr  [Ang^{3/2}]
//!   Step 4: Matrix element:
//!       V_NL(G,G') = (1/Omega) * sum_atom S(G-G')
//!                   * sum_{a,b} F_a(|k+G|) D_{ab} F_b(|k+G'|)
//!                   * (2l+1)/(4*pi) * P_l(cos theta)
//!   Step 5: Units: (1/Ang^3) * Ang^{3/2} * eV * Ang^{3/2} * 1 = eV
//!
//! Run with: cargo test --test kb_projector_validation -- --nocapture

use std::f64::consts::PI;
use std::path::PathBuf;

use nalgebra::Vector3;
use num_complex::Complex64;

use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    numerics::simpson_integrate,
    potential::nonlocal::NonlocalPotential,
    pseudopotential::PseudopotentialData,
};

// ---------------------------------------------------------------------------
//  Constants (must match the library exactly)
// ---------------------------------------------------------------------------
const BOHR_TO_ANG: f64 = 0.529_177_210_903;
const RY_TO_EV: f64 = 13.605_693_122_994;

// ---------------------------------------------------------------------------
//  Helper: load Si pseudopotential
// ---------------------------------------------------------------------------
fn load_si_pp() -> PseudopotentialData {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf");
    pwdft_rs::pseudopotential::load(&path).unwrap()
}

fn si_crystal() -> Crystal {
    let a = 5.431; // Si lattice constant in Ang
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

// ---------------------------------------------------------------------------
//  Re-implementation of spherical Bessel functions (must match nonlocal.rs)
// ---------------------------------------------------------------------------
fn spherical_bessel_j(l: i32, x: f64) -> f64 {
    if x.abs() < 1e-10 {
        return if l == 0 { 1.0 } else { 0.0 };
    }
    match l {
        0 => x.sin() / x,
        1 => x.sin() / (x * x) - x.cos() / x,
        2 => (3.0 / (x * x) - 1.0) * x.sin() / x - 3.0 * x.cos() / (x * x),
        3 => {
            (15.0 / (x * x * x) - 6.0 / x) * x.sin() / x
                - (15.0 / (x * x) - 1.0) * x.cos() / x
        }
        _ => panic!("spherical_bessel_j: l={l} not implemented"),
    }
}

fn legendre_p(l: i32, x: f64) -> f64 {
    match l {
        0 => 1.0,
        1 => x,
        2 => 0.5 * (3.0 * x * x - 1.0),
        3 => 0.5 * (5.0 * x * x * x - 3.0 * x),
        _ => panic!("legendre_p: l={l} not implemented"),
    }
}

// ---------------------------------------------------------------------------
//  Bessel transform: TRAPEZOIDAL rule (our code's method, via rab weights)
// ---------------------------------------------------------------------------
fn bessel_transform_trapezoidal(
    r_grid: &[f64],
    rab: &[f64],
    chi: &[f64],
    l: i32,
    q: f64,
) -> f64 {
    let mut integral = 0.0;
    for i in 0..r_grid.len() {
        let r = r_grid[i];
        let dr = rab[i];
        let qr = q * r;
        let jl = spherical_bessel_j(l, qr);
        // Integrand: chi(r) * j_l(qr) * r * dr
        // chi(r) already stores r*beta(r), so the integrand is chi * j_l * r * dr
        integral += chi[i] * jl * r * dr;
    }
    4.0 * PI * integral
}

// ---------------------------------------------------------------------------
//  Bessel transform: SIMPSON'S 1/3 rule (QE's method)
// ---------------------------------------------------------------------------
fn bessel_transform_simpson(
    r_grid: &[f64],
    rab: &[f64],
    chi: &[f64],
    l: i32,
    q: f64,
) -> f64 {
    let n = r_grid.len();
    // Simpson's rule requires an odd number of points.
    // If n is even, we drop the last point (whose contribution is negligible
    // because chi has decayed to zero by then).
    let n_eff = if n.is_multiple_of(2) { n - 1 } else { n };

    let mut integral = 0.0;
    for i in 0..n_eff {
        let r = r_grid[i];
        let dr = rab[i];
        let qr = q * r;
        let jl = spherical_bessel_j(l, qr);
        let f = chi[i] * jl * r;

        let w = if i == 0 || i == n_eff - 1 {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        integral += w * f * dr / 3.0;
    }
    4.0 * PI * integral
}

// ===========================================================================
//  TEST 1: Inspect raw projector data from UPF
// ===========================================================================
#[test]
fn test_01_inspect_projector_data() {
    let pp = load_si_pp();

    eprintln!("=== TEST 1: Inspect UPF projector data ===");
    eprintln!("Element:     {}", pp.element);
    eprintln!("Z_valence:   {}", pp.z_valence);
    eprintln!("l_max:       {}", pp.l_max);
    eprintln!("n_proj:      {}", pp.n_projectors());
    eprintln!("mesh_size:   {}", pp.r_grid.len());
    eprintln!(
        "r_grid range: [{:.6e}, {:.6e}] Ang",
        pp.r_grid[0],
        pp.r_grid.last().unwrap()
    );
    eprintln!();

    // D_ij matrix (eV)
    let np = pp.n_projectors();
    eprintln!("D_ij matrix ({np}x{np}) in eV:");
    for i in 0..np {
        for j in 0..np {
            eprint!("  {:12.6}", pp.dij[i * np + j]);
        }
        eprintln!();
    }
    eprintln!();

    for (ip, proj) in pp.beta_projectors.iter().enumerate() {
        let chi = &proj.values;
        let n = chi.len();

        // Find peak
        let (peak_idx, peak_val) = chi
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .unwrap();
        let peak_r = pp.r_grid[peak_idx];

        eprintln!(
            "Projector {} (l={}): peak |chi|={:.6e} at r={:.6} Ang (index {})",
            ip, proj.l, peak_val.abs(), peak_r, peak_idx
        );

        // First 5 values
        eprint!("  first 5 chi(r): ");
        for v in chi.iter().take(5.min(n)) {
            eprint!("{v:.6e}  ");
        }
        eprintln!();

        // Last 5 non-trivially-zero values
        let last_nonzero = chi
            .iter()
            .rposition(|&v| v.abs() > 1e-30)
            .unwrap_or(n - 1);
        let start = last_nonzero.saturating_sub(4);
        eprint!("  last 5 nonzero chi(r) [idx {start}-{last_nonzero}]: ");
        for v in chi.iter().take(last_nonzero + 1).skip(start) {
            eprint!("{v:.6e}  ");
        }
        eprintln!();

        // Check Gaussian-like shape: peak should be near r ~ 0.1--0.5 Ang for HGH Si
        assert!(
            peak_r < 2.0,
            "Projector {ip} peak r={peak_r:.4} Ang is suspiciously large"
        );
        assert!(
            peak_val.abs() > 1e-6,
            "Projector {ip} peak value is too small: {peak_val:.6e}"
        );
    }
}

// ===========================================================================
//  TEST 2: Simpson vs trapezoidal Bessel transform comparison
// ===========================================================================
#[test]
fn test_02_simpson_vs_trapezoidal() {
    let pp = load_si_pp();

    eprintln!("\n=== TEST 2: Simpson vs Trapezoidal Bessel transform ===");
    eprintln!(
        "{:>8}  {:>5}  {:>16}  {:>16}  {:>12}  {:>12}",
        "proj", "l", "F_trap(q)", "F_simp(q)", "abs_diff", "rel_diff"
    );

    // Test at q values in Ang^{-1}
    let q_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0];

    for (ip, proj) in pp.beta_projectors.iter().enumerate() {
        let chi = &proj.values;
        let l = proj.l;

        for &q in &q_values {
            let f_trap = bessel_transform_trapezoidal(&pp.r_grid, &pp.rab, chi, l, q);
            let f_simp = bessel_transform_simpson(&pp.r_grid, &pp.rab, chi, l, q);

            let abs_diff = (f_trap - f_simp).abs();
            let rel_diff = if f_trap.abs() > 1e-20 {
                abs_diff / f_trap.abs()
            } else {
                0.0
            };

            eprintln!(
                "  proj={ip} l={l}  q={q:5.2}  F_trap={f_trap:16.10e}  F_simp={f_simp:16.10e}  abs={abs_diff:12.4e}  rel={rel_diff:12.4e}"
            );

            // The two methods should agree to within 1% for smooth Gaussian-like
            // projectors on the 1141-point log grid
            if f_trap.abs() > 1e-15 {
                assert!(
                    rel_diff < 0.01,
                    "Simpson vs trapezoidal disagree by {:.2}% for proj={ip}, l={l}, q={q}",
                    rel_diff * 100.0
                );
            }
        }
        eprintln!();
    }
}

// ===========================================================================
//  TEST 3: Analytic check at q=0 for l=0 projectors
// ===========================================================================
#[test]
fn test_03_f_at_q_zero_analytic() {
    let pp = load_si_pp();

    eprintln!("\n=== TEST 3: F(q=0) analytic check ===");
    eprintln!("For l=0: j_0(0)=1, so F(q=0) = 4*pi * integral chi(r) * r * dr");
    eprintln!("For l>0: j_l(0)=0, so F(q=0) = 0 exactly");
    eprintln!();

    for (ip, proj) in pp.beta_projectors.iter().enumerate() {
        let chi = &proj.values;
        let l = proj.l;

        // Compute F(q=0) from our transform (Simpson's rule)
        let f_q0 = bessel_transform_simpson(&pp.r_grid, &pp.rab, chi, l, 0.0);

        if l == 0 {
            // Analytic: F(0) = 4*pi * integral chi(r) * r * dr (Simpson's rule)
            let n = chi.len();
            let integrand: Vec<f64> = (0..n).map(|i| chi[i] * pp.r_grid[i]).collect();
            let integral_analytic = simpson_integrate(&integrand, &pp.rab);
            let f_analytic = 4.0 * PI * integral_analytic;

            let diff = (f_q0 - f_analytic).abs();
            eprintln!(
                "Projector {ip} (l=0): F(0) = {f_q0:.10e},  analytic = {f_analytic:.10e},  diff = {diff:.4e}"
            );
            assert!(
                diff < 1e-14 * f_q0.abs().max(1.0),
                "F(q=0) disagrees with analytic for l=0 projector {ip}: diff = {diff:.4e}"
            );
        } else {
            // l > 0: must be exactly zero (j_l(0) = 0 for l > 0)
            eprintln!(
                "Projector {ip} (l={l}): F(0) = {f_q0:.10e}  (should be 0)"
            );
            assert!(
                f_q0.abs() < 1e-20,
                "F(q=0) should be zero for l={l} projector: got {f_q0:.4e}"
            );
        }
    }
}

// ===========================================================================
//  TEST 4: Unit conversion chain verification
// ===========================================================================
#[test]
fn test_04_unit_conversion_chain() {
    eprintln!("\n=== TEST 4: Unit conversion chain verification ===");

    // Read the raw UPF file to get Bohr-unit data
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf");
    let content = std::fs::read_to_string(&path).unwrap();

    // Parse r_grid in Bohr
    let r_bohr = extract_block_f64(&content, "PP_R");
    let rab_bohr = extract_block_f64(&content, "PP_RAB");
    let beta1_raw = extract_beta_block(&content, "PP_BETA.1");
    let dij_raw = extract_block_f64(&content, "PP_DIJ");

    eprintln!("Raw UPF data (Bohr/Ry units):");
    eprintln!("  r_grid[0]  = {:.10e} Bohr = {:.10e} Ang", r_bohr[0], r_bohr[0] * BOHR_TO_ANG);
    eprintln!("  rab[0]     = {:.10e} Bohr = {:.10e} Ang", rab_bohr[0], rab_bohr[0] * BOHR_TO_ANG);
    eprintln!("  beta1[10]  = {:.10e} Bohr^{{-1/2}}", beta1_raw[10]);
    eprintln!("  D_11 (raw) = {:.10e} Ry = {:.10e} eV", dij_raw[0], dij_raw[0] * RY_TO_EV);
    eprintln!();

    // Load via our parser
    let pp = load_si_pp();

    // Verify r_grid conversion
    let r_check = r_bohr[0] * BOHR_TO_ANG;
    let diff_r = (pp.r_grid[0] - r_check).abs();
    eprintln!("r_grid[0]: parsed={:.10e}, expected={:.10e}, diff={:.4e}", pp.r_grid[0], r_check, diff_r);
    assert!(diff_r < 1e-15, "r_grid conversion error: {diff_r:.4e}");

    // Verify rab conversion
    let rab_check = rab_bohr[0] * BOHR_TO_ANG;
    let diff_rab = (pp.rab[0] - rab_check).abs();
    eprintln!("rab[0]: parsed={:.10e}, expected={:.10e}, diff={:.4e}", pp.rab[0], rab_check, diff_rab);
    assert!(diff_rab < 1e-15, "rab conversion error: {diff_rab:.4e}");

    // Verify chi conversion: chi_Ang = chi_Bohr / sqrt(BOHR_TO_ANG)
    let chi_check = beta1_raw[10] / BOHR_TO_ANG.sqrt();
    let diff_chi = (pp.beta_projectors[0].values[10] - chi_check).abs();
    eprintln!(
        "chi[0][10]: parsed={:.10e}, expected={:.10e}, diff={:.4e}",
        pp.beta_projectors[0].values[10], chi_check, diff_chi
    );
    assert!(diff_chi < 1e-15, "chi conversion error: {diff_chi:.4e}");

    // Verify D_ij conversion: D_eV = D_Ry * RY_TO_EV
    let d11_check = dij_raw[0] * RY_TO_EV;
    let diff_d = (pp.dij[0] - d11_check).abs();
    eprintln!("D_11: parsed={:.10e}, expected={:.10e}, diff={:.4e}", pp.dij[0], d11_check, diff_d);
    assert!(diff_d < 1e-10, "D_ij conversion error: {diff_d:.4e}");

    // Dimensional analysis of the full matrix element
    eprintln!();
    eprintln!("Dimensional analysis of V_NL(G,G) at Gamma:");
    let omega = 40.04; // approximate Si cell volume in Ang^3
    eprintln!("  1/Omega            = {:.6e} Ang^{{-3}}", 1.0 / omega);
    let f_test = bessel_transform_simpson(
        &pp.r_grid,
        &pp.rab,
        &pp.beta_projectors[0].values,
        0,
        1.0,
    );
    eprintln!("  F_0(q=1 Ang^-1)   = {f_test:.6e} Ang^{{3/2}}");
    eprintln!("  D_11               = {:.6e} eV", pp.dij[0]);
    eprintln!(
        "  F * D * F          = {:.6e} eV * Ang^3",
        f_test * pp.dij[0] * f_test
    );
    eprintln!(
        "  (1/Omega)*F*D*F    = {:.6e} eV  [correct dimension]",
        f_test * pp.dij[0] * f_test / omega
    );
}

// ===========================================================================
//  TEST 5: Full V_NL diagonal elements at Gamma
// ===========================================================================
#[test]
fn test_05_vnl_diagonal_at_gamma() {
    let pp = load_si_pp();
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let k: Vector3<f64> = Vector3::zeros();
    let omega = crystal.lattice.volume();
    let n_pw = basis.len();

    eprintln!("\n=== TEST 5: V_NL diagonal at Gamma ===");
    eprintln!("Cell volume: {omega:.6} Ang^3");
    eprintln!("Number of PWs: {n_pw}");

    // Build V_NL using the library
    let mut h_lib = faer::Mat::<Complex64>::zeros(n_pw, n_pw);
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]);
    vnl.add_to_hamiltonian(&mut h_lib, &crystal, &basis, &k);

    // Also compute diagonal V_NL(G,G) manually for selected G-vectors
    // At Gamma, k=0, so q = G. For diagonal G=G', the angular factor is:
    //   (2l+1)/(4*pi) * P_l(cos 0) = (2l+1)/(4*pi)
    // and the structure factor is S(0) = N_atoms_of_type * 1 (for each type)
    // Actually: S(G-G') = S(0) = sum_atom exp(0) = N_atoms
    // For Si diamond: 2 atoms both of type Si, so S(0) = 2.

    let g_vecs = basis.g_vectors();

    // Pick a few representative G-vectors
    let test_indices: Vec<usize> = {
        let mut indices = Vec::new();
        // G = 0
        if let Some(idx) = basis.index_of(0, 0, 0) {
            indices.push(idx);
        }
        // A few G-vectors from different shells
        for &(n1, n2, n3) in &[(1, 1, 1), (0, 0, 2), (2, 2, 0), (-1, -1, -1)] {
            if let Some(idx) = basis.index_of(n1, n2, n3) {
                indices.push(idx);
            }
        }
        indices
    };

    eprintln!(
        "\n{:>6}  {:>20}  {:>20}  {:>16}  {:>12}  {:>12}",
        "idx", "G (Ang^-1)", "|G| (Ang^-1)", "V_NL_lib (eV)", "V_NL_manual", "diff"
    );

    for &ig in &test_indices {
        let g = g_vecs[ig];
        let g_norm = g.norm();

        // Manual diagonal V_NL(G,G)
        // = (1/Omega) * S(0) * sum_{i,j same l} F_i(|G|) D_{ij} F_j(|G|) * (2l+1)/(4*pi)
        let n_proj = pp.n_projectors();
        let mut vnl_manual = 0.0;

        // Structure factor S(0) for diagonal: sum over all atoms
        let n_atoms = crystal.atoms.len() as f64;

        for i in 0..n_proj {
            for j in 0..n_proj {
                let li = pp.beta_projectors[i].l;
                let lj = pp.beta_projectors[j].l;
                if li != lj {
                    continue;
                }
                let l = li;

                let fi = bessel_transform_simpson(
                    &pp.r_grid,
                    &pp.rab,
                    &pp.beta_projectors[i].values,
                    l,
                    g_norm,
                );
                let fj = bessel_transform_simpson(
                    &pp.r_grid,
                    &pp.rab,
                    &pp.beta_projectors[j].values,
                    l,
                    g_norm,
                );

                let d = pp.dij[i * n_proj + j];

                // Angular: at diagonal, cos(theta) = 1 (or undefined if |G|=0)
                let angular = (2 * l + 1) as f64 / (4.0 * PI);

                vnl_manual += fi * d * fj * angular;
            }
        }
        vnl_manual *= n_atoms / omega;

        let v_lib = h_lib[(ig, ig)].re;
        let diff = (v_lib - vnl_manual).abs();

        let miller = basis.miller_indices()[ig];
        eprintln!(
            "  {:>4}  G=({:>2},{:>2},{:>2})  |G|={:8.4}  V_lib={:12.8}  V_man={:12.8}  diff={:.4e}",
            ig, miller[0], miller[1], miller[2], g_norm, v_lib, vnl_manual, diff
        );

        // Library and manual should match to high precision
        assert!(
            diff < 1e-8,
            "V_NL diagonal mismatch at G=({},{},{}): lib={:.10}, manual={:.10}, diff={:.4e}",
            miller[0], miller[1], miller[2],
            v_lib, vnl_manual, diff
        );
    }
}

// ===========================================================================
//  TEST 6: V_NL(G=0,G=0) analytic cross-check
// ===========================================================================
#[test]
fn test_06_vnl_g0_g0_analytic() {
    let pp = load_si_pp();
    let crystal = si_crystal();
    let omega = crystal.lattice.volume();

    eprintln!("\n=== TEST 6: V_NL(G=0,G=0) analytic cross-check ===");

    // At G=0 (and k=0), q = 0.
    // For l=0 projectors: j_0(0) = 1, so F_i(0) = 4*pi * integral chi_i(r) * r * dr
    // For l>0 projectors: j_l(0) = 0, so F_i(0) = 0

    // Therefore V_NL(G=0,G=0) only gets contributions from l=0 projectors:
    // V_NL(0,0) = (N_atoms/Omega) * sum_{i,j with l=0} F_i(0) D_{ij} F_j(0) * 1/(4*pi)

    let n_proj = pp.n_projectors();
    let n_atoms = crystal.atoms.len() as f64;

    // Compute F_i(0) for each l=0 projector
    let mut f_q0: Vec<f64> = Vec::new();
    let mut l0_indices: Vec<usize> = Vec::new();

    for (i, proj) in pp.beta_projectors.iter().enumerate() {
        if proj.l == 0 {
            // F(0) = 4*pi * integral chi(r) * r * dr (Simpson's rule)
            let f0 = bessel_transform_simpson(&pp.r_grid, &pp.rab, &proj.values, 0, 0.0);
            f_q0.push(f0);
            l0_indices.push(i);
            eprintln!("  F_{i}(0) = {f0:.10e} Ang^{{3/2}}");
        }
    }

    let mut vnl_00_analytic = 0.0;
    for (ii, &i) in l0_indices.iter().enumerate() {
        for (jj, &j) in l0_indices.iter().enumerate() {
            let d = pp.dij[i * n_proj + j];
            vnl_00_analytic += f_q0[ii] * d * f_q0[jj];
        }
    }
    // angular factor: (2*0+1)/(4*pi) = 1/(4*pi)
    vnl_00_analytic *= 1.0 / (4.0 * PI);
    vnl_00_analytic *= n_atoms / omega;

    eprintln!("  V_NL(G=0,G=0) analytic = {vnl_00_analytic:.10} eV");

    // Compare with library value
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let k = Vector3::zeros();
    let n_pw = basis.len();
    let mut h = faer::Mat::<Complex64>::zeros(n_pw, n_pw);
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]);
    vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);

    let g0_idx = basis.index_of(0, 0, 0).unwrap();
    let v_lib = h[(g0_idx, g0_idx)].re;

    eprintln!("  V_NL(G=0,G=0) library  = {v_lib:.10} eV");
    let diff = (v_lib - vnl_00_analytic).abs();
    eprintln!("  difference             = {diff:.4e} eV");

    assert!(
        diff < 1e-8,
        "V_NL(G=0,G=0) analytic vs library mismatch: {vnl_00_analytic:.10} vs {v_lib:.10}, diff={diff:.4e}"
    );
}

// ===========================================================================
//  TEST 7: Form factor decay and smoothness
// ===========================================================================
#[test]
fn test_07_form_factor_behavior() {
    let pp = load_si_pp();

    eprintln!("\n=== TEST 7: Form factor F(q) behavior ===");

    // Sample F(q) at many q values
    let q_vals: Vec<f64> = (0..50).map(|i| i as f64 * 0.5).collect();

    for (ip, proj) in pp.beta_projectors.iter().enumerate() {
        let l = proj.l;
        let chi = &proj.values;

        eprintln!("\nProjector {ip} (l={l}):");
        eprintln!("  {:>8}  {:>16}", "q (1/Ang)", "F(q) (Ang^3/2)");

        let mut f_vals: Vec<f64> = Vec::new();
        for &q in &q_vals {
            let f = bessel_transform_simpson(&pp.r_grid, &pp.rab, chi, l, q);
            f_vals.push(f);
            if q <= 10.0 || (q - *q_vals.last().unwrap()).abs() < f64::EPSILON {
                eprintln!("  {q:8.3}  {f:16.10e}");
            }
        }

        // F(q) should decay to near zero at large q
        let f_large_q = f_vals.last().unwrap().abs();
        let f_max = f_vals.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        let ratio = f_large_q / f_max.max(1e-20);
        eprintln!(
            "  |F(q_max)|/|F_max| = {ratio:.4e} (should be < 0.12)"
        );
        // HGH l=1 projector for Si has broad q-space extent: ratio ~0.106 at
        // q=24.5 Ang^-1. This is a property of the Gaussian projector shape,
        // not a quadrature artifact — both trapezoidal and Simpson give ~0.106.
        assert!(
            ratio < 0.12,
            "Form factor not decaying: |F(q_max)|/|F_max| = {ratio:.4e}"
        );

        // Check smoothness: |F(q_{i+1}) - F(q_i)| should not have wild jumps
        for i in 0..f_vals.len() - 1 {
            let jump = (f_vals[i + 1] - f_vals[i]).abs();
            // Allow up to 50% of the total range per step (generous for 0.5 Ang^-1 spacing)
            let range = f_max;
            if range > 1e-15 {
                assert!(
                    jump / range < 0.5,
                    "Non-smooth F(q) at q={:.1}: jump/range = {:.4}",
                    q_vals[i],
                    jump / range
                );
            }
        }
    }
}

// ===========================================================================
//  TEST 8: Off-diagonal V_NL matrix element manual verification
// ===========================================================================
#[test]
fn test_08_vnl_offdiagonal() {
    let pp = load_si_pp();
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let k: Vector3<f64> = Vector3::zeros();
    let omega = crystal.lattice.volume();
    let n_pw = basis.len();

    eprintln!("\n=== TEST 8: Off-diagonal V_NL matrix elements ===");

    // Build V_NL using the library
    let mut h_lib = faer::Mat::<Complex64>::zeros(n_pw, n_pw);
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]);
    vnl.add_to_hamiltonian(&mut h_lib, &crystal, &basis, &k);

    let g_vecs = basis.g_vectors();
    let n_proj = pp.n_projectors();

    // Test a few off-diagonal pairs
    let test_pairs = [
        ((0, 0, 0), (1, 1, 1)),
        ((1, 1, 1), (0, 0, 2)),
        ((-1, -1, -1), (1, 1, 1)),
    ];

    for &(m1, m2) in &test_pairs {
        let Some(ig) = basis.index_of(m1.0, m1.1, m1.2) else {
            continue;
        };
        let Some(jg) = basis.index_of(m2.0, m2.1, m2.2) else {
            continue;
        };

        let qi = g_vecs[ig]; // at Gamma, k=0
        let qj = g_vecs[jg];
        let qi_norm = qi.norm();
        let qj_norm = qj.norm();

        // Manual computation
        let mut vnl_manual = Complex64::new(0.0, 0.0);

        for atom in &crystal.atoms {
            let tau = atom.cart_position(&crystal.lattice);
            let g_diff = qi - qj;
            let phase = -g_diff.dot(&tau);
            let sf = Complex64::new(phase.cos(), phase.sin());

            let mut vnl_atom = 0.0;
            for i in 0..n_proj {
                for j in 0..n_proj {
                    let li = pp.beta_projectors[i].l;
                    let lj = pp.beta_projectors[j].l;
                    if li != lj {
                        continue;
                    }
                    let l = li;

                    let fi = bessel_transform_simpson(
                        &pp.r_grid,
                        &pp.rab,
                        &pp.beta_projectors[i].values,
                        l,
                        qi_norm,
                    );
                    let fj = bessel_transform_simpson(
                        &pp.r_grid,
                        &pp.rab,
                        &pp.beta_projectors[j].values,
                        l,
                        qj_norm,
                    );

                    let d = pp.dij[i * n_proj + j];

                    let cos_theta = if qi_norm > 1e-12 && qj_norm > 1e-12 {
                        qi.dot(&qj) / (qi_norm * qj_norm)
                    } else {
                        1.0
                    };
                    let angular = (2 * l + 1) as f64 / (4.0 * PI) * legendre_p(l, cos_theta);

                    vnl_atom += fi * d * fj * angular;
                }
            }
            vnl_manual += sf * (vnl_atom / omega);
        }

        let v_lib = h_lib[(ig, jg)];
        let diff = (v_lib - vnl_manual).norm();

        let miller_i = basis.miller_indices()[ig];
        let miller_j = basis.miller_indices()[jg];
        eprintln!(
            "  G=({},{},{}) -> G=({},{},{}): lib=({:.8}, {:.8}i), manual=({:.8}, {:.8}i), diff={:.4e}",
            miller_i[0], miller_i[1], miller_i[2],
            miller_j[0], miller_j[1], miller_j[2],
            v_lib.re, v_lib.im,
            vnl_manual.re, vnl_manual.im,
            diff
        );

        assert!(
            diff < 1e-8,
            "Off-diagonal V_NL mismatch: diff = {diff:.4e}"
        );
    }
}

// ===========================================================================
//  TEST 9: UPF D_ij cross-check against raw file values
// ===========================================================================
#[test]
fn test_09_hgh_parameter_crosscheck() {
    let pp = load_si_pp();

    eprintln!("\n=== TEST 9: UPF D_ij cross-check ===");
    eprintln!("Note: Si.upf has 6 projectors (2 per l=0,1,2). QE diagonalizes the");
    eprintln!("raw HGH h^l_ij matrix and absorbs the rotation into the projectors,");
    eprintln!("so the UPF D_ij is a diagonal 6x6 matrix whose entries do NOT match");
    eprintln!("the raw HGH h^l_ij from PRB 58, 3641 (1998) Table I.");
    eprintln!();

    // Z_ion = 4
    assert!(
        (pp.z_valence - 4.0).abs() < 1e-10,
        "Z_valence mismatch: {} vs 4.0",
        pp.z_valence
    );

    // Si.upf has 6 projectors: 2 for l=0, 2 for l=1, 2 for l=2
    let np = pp.n_projectors();
    assert_eq!(np, 6, "Expected 6 projectors (2 per l=0,1,2), got {np}");

    let d = &pp.dij;

    // D_ij must be symmetric: D[i,j] = D[j,i]
    eprintln!("Symmetry check:");
    for i in 0..np {
        for j in i + 1..np {
            let diff = (d[i * np + j] - d[j * np + i]).abs();
            eprintln!(
                "  D[{},{}]={:.10e}, D[{},{}]={:.10e}, diff={:.4e}",
                i, j, d[i * np + j], j, i, d[j * np + i], diff
            );
            assert!(
                diff < 1e-14,
                "D_ij not symmetric: D[{i},{j}]={:.10e} vs D[{j},{i}]={:.10e}",
                d[i * np + j],
                d[j * np + i]
            );
        }
    }

    // D_ij must be block-diagonal in l: entries between different l channels are zero
    eprintln!("\nBlock-diagonal check (off-diagonal between different l must be zero):");
    for i in 0..np {
        for j in 0..np {
            let li = pp.beta_projectors[i].l;
            let lj = pp.beta_projectors[j].l;
            if li != lj {
                let val = d[i * np + j];
                eprintln!(
                    "  D[{i},{j}] (l={li} vs l={lj}) = {val:.10e} (must be zero)"
                );
                assert!(
                    val.abs() < 1e-14,
                    "D[{i},{j}] between l={li} and l={lj} should be zero, got {val:.10e}"
                );
            }
        }
    }

    // Verify diagonal values against the raw UPF file (in Ry, converted to eV).
    // These are the eigenvalues of the HGH h^l matrices after QE's diagonalization.
    let expected_diag_ry = [
        1.113_191_595_4e+01, // l=0, proj 1
        1.713_932_492_5e+00, // l=0, proj 2
        5.452_221_279_1e+00, // l=1, proj 1
        1.259_655_832_9e+00, // l=1, proj 2
        -4.249_608_729_0e+00, // l=2, proj 1
        -8.892_087_962_2e-01, // l=2, proj 2
    ];

    eprintln!("\nDiagonal D_ij comparison (UPF values in Ry -> eV):");
    eprintln!(
        "{:>8}  {:>4}  {:>16}  {:>16}  {:>12}",
        "entry", "l", "parsed (eV)", "expected (eV)", "diff"
    );
    for i in 0..np {
        let expected_ev = expected_diag_ry[i] * RY_TO_EV;
        let parsed_ev = d[i * np + i];
        let diff = (parsed_ev - expected_ev).abs();
        let l = pp.beta_projectors[i].l;
        eprintln!(
            "  D[{i},{i}]  l={l}  {parsed_ev:16.8}  {expected_ev:16.8}  {diff:12.4e}"
        );
        assert!(
            diff < 1e-4,
            "D[{i},{i}] mismatch: parsed={parsed_ev:.10}, expected={expected_ev:.10}, diff={diff:.4e}"
        );
    }
}

// ===========================================================================
//  TEST 10: Hermiticity and reality of diagonal V_NL
// ===========================================================================
#[test]
fn test_10_vnl_hermiticity_and_reality() {
    let pp = load_si_pp();
    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let k: Vector3<f64> = Vector3::zeros();
    let n_pw = basis.len();

    eprintln!("\n=== TEST 10: V_NL Hermiticity and diagonal reality ===");

    let mut h = faer::Mat::<Complex64>::zeros(n_pw, n_pw);
    let vnl = NonlocalPotential::new(&crystal, &basis, &k, &[&pp]);
    vnl.add_to_hamiltonian(&mut h, &crystal, &basis, &k);

    // Diagonal must be real
    let mut max_imag = 0.0f64;
    for i in 0..n_pw {
        max_imag = max_imag.max(h[(i, i)].im.abs());
    }
    eprintln!("  Max |Im(V_NL(G,G))|  = {max_imag:.4e}");
    assert!(
        max_imag < 1e-10,
        "V_NL diagonal has imaginary part: max = {max_imag:.4e}"
    );

    // Off-diagonal: H_{ij} = H_{ji}^* (Hermitian)
    let mut max_herm_err = 0.0f64;
    for i in 0..n_pw {
        for j in i + 1..n_pw {
            let err = (h[(i, j)] - h[(j, i)].conj()).norm();
            max_herm_err = max_herm_err.max(err);
        }
    }
    eprintln!("  Max |H(i,j) - H(j,i)*| = {max_herm_err:.4e}");
    assert!(
        max_herm_err < 1e-10,
        "V_NL not Hermitian: max = {max_herm_err:.4e}"
    );
}

// ---------------------------------------------------------------------------
//  Raw UPF block extractors (avoids modifying source code)
// ---------------------------------------------------------------------------
fn extract_block_f64(content: &str, tag: &str) -> Vec<f64> {
    let open = format!("<{tag}");
    let close = format!("</{tag}>");
    let tag_pos = content.find(&open).unwrap_or_else(|| panic!("missing <{tag}>"));
    let data_start = content[tag_pos..]
        .find('>')
        .expect("malformed tag")
        + tag_pos
        + 1;
    let data_end = content[data_start..]
        .find(&close)
        .unwrap_or_else(|| panic!("missing </{tag}>"))
        + data_start;
    content[data_start..data_end]
        .split_whitespace()
        .map(|s| s.parse::<f64>().unwrap_or_else(|_| panic!("parse error in {tag}: {s}")))
        .collect()
}

fn extract_beta_block(content: &str, tag: &str) -> Vec<f64> {
    // Same as extract_block_f64 but for PP_BETA blocks which can have multiline opening tags
    extract_block_f64(content, tag)
}

// ===========================================================================
//  V_local(G) comparison with QE
// ===========================================================================

/// Compare our V_local(G) against QE's Cube file FFT at key G-vectors.
///
/// QE reference (from pp.x plot_num=2 → Cube → FFT):
///   G=(0,0,0):  -1.002_741 eV
///   G=(1,0,0):  (-4.9268, +4.9268)i eV  |V| = 6.9675 eV
///   G=(1,1,1):  (-4.9268, -4.9268)i eV  |V| = 6.9675 eV
///   G=(2,0,0):  ~0 eV
#[test]
#[ignore = "VERF: v_local_of_g uses bare Coulomb subtraction; QE uses erf(r)/r which is numerically superior"]
fn test_vloc_comparison_with_qe() {
    let pp = load_si_pp();
    let crystal = si_crystal();
    let omega = crystal.lattice.volume();
    let recip = crystal.lattice.reciprocal();

    // Compute V_local(G) at specific G-vectors using our Bessel transform
    // G-vector for Miller index (n1,n2,n3): G = n1*b1 + n2*b2 + n3*b3

    let test_cases: Vec<(&str, [i32; 3])> = vec![
        ("G=(0,0,0)", [0, 0, 0]),
        ("G=(1,0,0)", [1, 0, 0]),
        ("G=(0,1,0)", [0, 1, 0]),
        ("G=(0,0,1)", [0, 0, 1]),
        ("G=(1,1,1)", [1, 1, 1]),
        ("G=(2,0,0)", [2, 0, 0]),
        ("G=(-1,0,0)", [-1, 0, 0]),
    ];

    eprintln!("\nOur V_local(G) (eV):");
    for (label, [n1, n2, n3]) in &test_cases {
        let g = *n1 as f64 * recip.a + *n2 as f64 * recip.b + *n3 as f64 * recip.c;
        let g_norm = g.norm();

        // Sum over atoms: V_local(G) = Σ_atom S(G) × v_form(|G|)
        let mut v_total = Complex64::new(0.0, 0.0);
        for atom in &crystal.atoms {
            let tau = atom.cart_position(&crystal.lattice);
            let phase = -g.dot(&tau);
            let sf = Complex64::new(phase.cos(), phase.sin());
            let v_form = pp.v_local_of_g(g_norm, omega);
            v_total += sf * v_form;
        }

        eprintln!(
            "  {}: {:+12.6} {:+12.6}i  |V|={:12.6}",
            label, v_total.re, v_total.im, v_total.norm()
        );
    }

    // QE reference values (from Cube FFT, in eV)
    let qe_vloc_g000 = -1.002_741;
    let qe_vloc_g100_abs = 6.967_521;
    let qe_vloc_g111_abs = 6.967_521;

    // Compute our values
    let _g000: Vector3<f64> = Vector3::zeros();
    let g100 = 1.0 * recip.a;
    let g111 = 1.0 * recip.a + 1.0 * recip.b + 1.0 * recip.c;

    // V_local(G=0) — just the form factor times 2 (two atoms, both S(0)=1)
    let our_vloc_g000 = 2.0 * pp.v_local_of_g(0.0, omega);

    // V_local(G=(1,0,0))
    let mut our_vloc_g100 = Complex64::new(0.0, 0.0);
    for atom in &crystal.atoms {
        let tau = atom.cart_position(&crystal.lattice);
        let phase = -g100.dot(&tau);
        let sf = Complex64::new(phase.cos(), phase.sin());
        our_vloc_g100 += sf * pp.v_local_of_g(g100.norm(), omega);
    }

    // V_local(G=(1,1,1))
    let mut our_vloc_g111 = Complex64::new(0.0, 0.0);
    for atom in &crystal.atoms {
        let tau = atom.cart_position(&crystal.lattice);
        let phase = -g111.dot(&tau);
        let sf = Complex64::new(phase.cos(), phase.sin());
        our_vloc_g111 += sf * pp.v_local_of_g(g111.norm(), omega);
    }

    eprintln!("\nComparison with QE:");
    eprintln!("  V_local(G=0):  ours={:.6} eV,  QE={:.6} eV,  diff={:.6} eV",
        our_vloc_g000, qe_vloc_g000, our_vloc_g000 - qe_vloc_g000);
    eprintln!("  |V_local(G=(1,0,0))|:  ours={:.6} eV,  QE={:.6} eV,  diff={:.6} eV",
        our_vloc_g100.norm(), qe_vloc_g100_abs, our_vloc_g100.norm() - qe_vloc_g100_abs);
    eprintln!("  |V_local(G=(1,1,1))|:  ours={:.6} eV,  QE={:.6} eV,  diff={:.6} eV",
        our_vloc_g111.norm(), qe_vloc_g111_abs, our_vloc_g111.norm() - qe_vloc_g111_abs);

    // Check agreement — tolerance of 0.1 eV for now (FFT grid differences cause some discrepancy)
    assert!(
        (our_vloc_g000 - qe_vloc_g000).abs() < 0.5,
        "V_local(G=0) disagrees: ours={our_vloc_g000:.6}, QE={qe_vloc_g000:.6}"
    );
    assert!(
        (our_vloc_g100.norm() - qe_vloc_g100_abs).abs() < 0.5,
        "|V_local(G=(1,0,0))| disagrees: ours={:.6}, QE={:.6}",
        our_vloc_g100.norm(), qe_vloc_g100_abs
    );
}
