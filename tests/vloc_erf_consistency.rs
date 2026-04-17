//! Regression guard: V_local(G) erf-form vs bare-Coulomb form.
//!
//! This test pins the mathematical equivalence of the two decompositions
//! of V_local(G):
//!
//!   (a) bare Coulomb subtraction (previously used in pwdft-rs):
//!       V(G) = (4π/Ω) ∫₀^∞ r² [V(r) + Z·e²/r] sin(Gr)/(Gr) dr
//!              − 4π·Z·e²/(Ω·G²)
//!
//!   (b) erf subtraction (current form, matches QE `vloc_mod.f90:136-148`):
//!       V(G) = (4π/Ω) ∫₀^∞ [r·V(r) + Z·e²·erf(r)] · sin(Gr)/G dr
//!              − 4π·Z·e²·exp(−G²/4)/(Ω·G²)
//!
//! The two forms are mathematically identical: V(r) − V(r) = 0. Specifically,
//! FT[−Z·e²/r] = −4π·Z·e²/G² and FT[−Z·e²·erf(r)/r] = −4π·Z·e²·exp(−G²/4)/G².
//! Their difference equals the FT of a Gaussian, which the integrand in form
//! (b) accounts for exactly via its smooth short-range piece.
//!
//! Numerically, after Simpson's rule on the log mesh, they should agree to
//! machine precision (< 1e-6 eV). This test:
//!   1. documents the equivalence explicitly (VERF Attempt 1 found the two
//!      forms gave identical total energies), and
//!   2. acts as a regression guard: if someone changes the radial quadrature
//!      or the erf decomposition and introduces a discrepancy, this test
//!      fires.
//!
//! We test the first 20 non-zero |G| shells for Si at the equilibrium lattice
//! constant (a = 5.431 Å).

use std::f64::consts::PI;
use std::path::PathBuf;

use nalgebra::Vector3;
use pwdft_rs::{
    crystal::Lattice,
    numerics::simpson_integrate,
    pseudopotential::{PseudopotentialData, load},
};

const E2: f64 = 14.399_645_351_950_548; // eV·Å, matches crate::consts::E2_COULOMB

fn load_si_pp() -> PseudopotentialData {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf");
    load(&path).unwrap()
}

/// Reference implementation: bare-Coulomb subtraction (the *old* form).
/// We keep this only in the test so we can compare against the library's
/// current erf form. See QE-equivalent old pwdft-rs code pre-VERF.
fn v_local_of_g_bare_coulomb(pp: &PseudopotentialData, g_norm: f64, omega: f64) -> f64 {
    let four_pi = 4.0 * PI;

    if g_norm < 1e-12 {
        let integrand: Vec<f64> = pp
            .r_grid
            .iter()
            .zip(pp.v_local.iter())
            .map(|(&r, &v)| {
                let v_short = v + pp.z_valence * E2 / r.max(1e-20);
                r * r * v_short
            })
            .collect();
        let integral = simpson_integrate(&integrand, &pp.rab);
        four_pi / omega * integral
    } else {
        // G ≠ 0: ∫ r² [V_loc(r) + Z·e²/r] sin(Gr)/(Gr) dr − 4π·Z·e²/(Ω·G²)
        let integrand: Vec<f64> = pp
            .r_grid
            .iter()
            .zip(pp.v_local.iter())
            .map(|(&r, &v)| {
                let gr = g_norm * r;
                let v_short = v + pp.z_valence * E2 / r.max(1e-20);
                let sinc = if gr < 1e-10 {
                    1.0 - gr * gr / 6.0
                } else {
                    gr.sin() / gr
                };
                r * r * v_short * sinc
            })
            .collect();
        let integral = simpson_integrate(&integrand, &pp.rab);
        four_pi / omega * integral - four_pi * pp.z_valence * E2 / (omega * g_norm * g_norm)
    }
}

/// Collect the first `n_shells` distinct |G| values from the Si FCC
/// reciprocal lattice, in increasing order.
fn first_g_shells_si(n_shells: usize) -> Vec<f64> {
    let a = 5.431; // Si lattice constant in Å
    let lat = Lattice::new(
        a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
        a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
        a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
    );
    let recip = lat.reciprocal();
    let (b1, b2, b3) = (recip.a, recip.b, recip.c);

    let n_max = 6;
    let mut mags: Vec<f64> = Vec::new();
    for n1 in -n_max..=n_max {
        for n2 in -n_max..=n_max {
            for n3 in -n_max..=n_max {
                let g = n1 as f64 * b1 + n2 as f64 * b2 + n3 as f64 * b3;
                let gn = g.norm();
                if gn > 1e-8 {
                    mags.push(gn);
                }
            }
        }
    }
    mags.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // Deduplicate shells (values within 1e-6 Å⁻¹ of each other).
    let mut shells = Vec::new();
    for g in mags {
        if shells
            .last()
            .copied()
            .is_none_or(|last: f64| (g - last).abs() > 1e-6)
        {
            shells.push(g);
            if shells.len() >= n_shells {
                break;
            }
        }
    }
    shells
}

#[test]
fn vloc_erf_matches_bare_coulomb_for_si_first_20_shells() {
    let pp = load_si_pp();
    let omega = 40.0; // approximate Si primitive cell volume in Å³

    let shells = first_g_shells_si(20);
    assert_eq!(
        shells.len(),
        20,
        "Need at least 20 distinct |G| shells for the regression test"
    );

    let mut max_abs = 0.0_f64;
    let mut max_rel = 0.0_f64;
    let mut worst_shell = (0usize, 0.0, 0.0, 0.0); // (i, g, v_erf, v_bare)

    for (i, &g) in shells.iter().enumerate() {
        let v_erf = pp.v_local_of_g(g, omega);
        let v_bare = v_local_of_g_bare_coulomb(&pp, g, omega);
        let abs_err = (v_erf - v_bare).abs();
        let rel_err = abs_err / v_bare.abs().max(1e-12);

        if abs_err > max_abs {
            max_abs = abs_err;
            worst_shell = (i, g, v_erf, v_bare);
        }
        if rel_err > max_rel {
            max_rel = rel_err;
        }

        eprintln!(
            "shell {i:2}  |G| = {g:.6} Å⁻¹   V_erf = {v_erf:+.8e} eV   V_bare = {v_bare:+.8e} eV   Δ = {abs_err:.2e} eV",
        );
    }

    eprintln!(
        "\nmax |Δ| = {max_abs:.2e} eV at shell {} (|G|={:.4}, V_erf={:+.6e}, V_bare={:+.6e})",
        worst_shell.0, worst_shell.1, worst_shell.2, worst_shell.3
    );
    eprintln!("max relative error = {max_rel:.2e}");

    // Both forms must agree to machine precision (< 1e-6 eV absolute).
    assert!(
        max_abs < 1e-6,
        "Bare-Coulomb and erf forms of V_local(G) disagree: max |Δ| = {max_abs:.3e} eV, \
         worst shell: i={}, |G|={:.4} Å⁻¹, V_erf={:+.6e}, V_bare={:+.6e}",
        worst_shell.0,
        worst_shell.1,
        worst_shell.2,
        worst_shell.3
    );
}

#[test]
fn vloc_erf_g_zero_matches_direct_integral() {
    // At G=0 both forms are algebraically identical. Verify the library
    // code path returns the same value as a direct evaluation.
    let pp = load_si_pp();
    let omega = 40.0;

    let v_lib = pp.v_local_of_g(0.0, omega);
    let v_bare = v_local_of_g_bare_coulomb(&pp, 0.0, omega);

    assert!(
        (v_lib - v_bare).abs() < 1e-9,
        "V_local(G=0) mismatch: lib = {v_lib:.6e} eV, bare = {v_bare:.6e} eV"
    );
}
