//! Parser for UPF v2 pseudopotential files (Quantum ESPRESSO format).
//!
//! UPF files use Rydberg atomic units: energies in Ry, lengths in Bohr.
//! We convert to internal units (eV, Å) on parse.
//!
//! Layout:
//! - `xml` — text-level helpers that pull attribute values and numeric data blocks out of the UPF
//!   XML.
//! - `convert` — unit-conversion body that assembles a
//!   [`crate::pseudopotential::PseudopotentialData`] from the parsed text.
//!
//! Only [`parse`] is public outside this folder; the helpers are
//! `pub(super)` and must not leak.

mod parse;
mod xml;

use std::path::Path;

use elements_rs::Element;
use parse::parse_upf_body;

use crate::{
    consts::{E2_COULOMB as E2, G_ZERO_THRESHOLD},
    error::{PwdftError, Result},
    numerics::simpson_integrate,
};

/// Unit-converted pseudopotential data in internal units (eV, Å).
///
/// Parsed from UPF v2 format.
#[derive(Debug, Clone)]
pub struct UpfPseudoPotential {
    /// Element symbol (e.g. "Si").
    pub element: Element,
    /// Number of valence electrons.
    pub z_valence: f64,
    /// Maximum angular momentum in non-local projectors.
    pub l_max: i32,
    /// Radial grid points (Å).
    pub r_grid: Vec<f64>,
    /// Integration weights dr for Simpson/trapezoidal rule (Å).
    pub rab: Vec<f64>,
    /// Local pseudopotential on radial grid (eV).
    /// Includes the -Z_val e²/r Coulomb tail.
    pub v_local: Vec<f64>,
    /// Non-local projectors: `beta[proj_index]` = (angular_momentum,
    /// radial_values). Radial values store χ(r) = r·β(r) in Å^{-1/2} (no
    /// energy dimension). Energy enters through D_ij. The KB matrix element
    /// is: V_NL = (1/Ω) Σ F_i D_ij F_j × angular, where F = 4π ∫ χ(r)
    /// j_l(qr) r dr.
    pub beta_projectors: Vec<BetaProjector>,
    /// D_ij coupling matrix for non-local projectors (eV).
    /// Stored as a flat n_proj × n_proj matrix in row-major order.
    pub dij: Vec<f64>,
    /// Atomic charge density on radial grid (e/Å, stores 4πr²ρ(r)).
    /// May be empty if not provided by the pseudopotential.
    pub rho_atom: Vec<f64>,
    /// Nonlinear core correction (NLCC) charge density on radial grid.
    ///
    /// Stores the **bare volumetric** density ρ_core(r) in e/Å³ (NOT
    /// 4πr²·ρ, which is the `PP_RHOATOM` convention — see the
    /// `rho_atom` field above). The downstream radial Bessel transform
    /// in `crate::scf::potentials::compute_core_density` multiplies by
    /// r² and 4π.
    ///
    /// Empty if the pseudopotential has no NLCC (`core_correction="F"`).
    ///
    /// Physics: ρ_core enters only `E_xc[ρ_val + ρ_core]` and
    /// `v_xc[ρ_val + ρ_core]`. It is *not* added to the Hartree source,
    /// is *not* counted as valence (no contribution to the electron
    /// count), and in LSDA is split symmetrically as ρ_core/2 per spin
    /// channel. See `crate::scf::energy::add_core_density` and
    /// Louie, Froyen, Cohen, *Phys. Rev. B* **26**, 1738 (1982).
    pub core_charge: Vec<f64>,
}

/// A single non-local beta projector.
#[derive(Debug, Clone)]
pub struct BetaProjector {
    /// Angular momentum quantum number.
    pub l: i32,
    /// Radial function β(r) on the radial grid.
    pub values: Vec<f64>,
}

impl UpfPseudoPotential {
    /// Load and parse a UPF v2 pseudopotential file.
    ///
    /// # Errors
    ///
    /// - `PwdftError::Io` (from the `?` on `std::fs::read_to_string`) if the file cannot be opened
    ///   or read.
    /// - `PwdftError::Parse` if the file does not look like UPF (no `<UPF` or `<PP_HEADER` marker)
    ///   — the loader rejects unknown formats rather than guessing.
    /// - Any `PwdftError::Parse` forwarded from [`upf::parse`] when the UPF body is malformed; see
    ///   that function's `# Errors` for the full list.
    pub fn load_from_path(path: &Path) -> Result<UpfPseudoPotential> {
        let content = std::fs::read_to_string(path)?;

        if content.contains("<UPF") || content.contains("<PP_HEADER") {
            parse_upf_body(&content)
        } else {
            Err(PwdftError::Parse(format!(
                "unrecognized pseudopotential format in {} (only UPF v2 is supported)",
                path.display()
            )))
        }
    }

    /// Load a UPF v2 pseudopotential from the workspace's library.
    ///
    /// Accepts an `Element`, a `u8` (atomic number), or a `&str` (symbol).
    pub fn load<E>(element_like: E) -> Result<Self>
    where
        E: TryInto<Element>,
        <E as TryInto<Element>>::Error: std::fmt::Display,
    {
        let element: Element = element_like
            .try_into()
            .map_err(|e| PwdftError::Parse(format!("Invalid element for pseudopotential: {e}")))?;

        let path = super::PP_LIBRARY_ROOT
            .join("nc/lda")
            .join(format!("{}.upf", element.symbol()));

        Self::load_from_path(&path)
    }

    /// Number of non-local beta projectors.
    pub fn n_projectors(&self) -> usize {
        self.beta_projectors.len()
    }

    /// Whether this PP has nonlinear core correction.
    pub fn has_nlcc(&self) -> bool {
        !self.core_charge.is_empty()
    }

    /// Compute V_local(G) via spherical Bessel transform.
    ///
    /// For G = 0:
    ///   V_local(0) = (4π/Ω) ∫₀^∞ r² [V_local(r) + Z e²/r] dr
    ///
    /// For G ≠ 0 we use the erf-subtracted form:
    ///   V_local(G) = (4π/Ω) ∫₀^∞ [r·V_local(r) + Z e²·erf(r)] · sin(Gr)/G dr
    ///                - (4π/Ω) Z e² · exp(-G²/4) / G²
    ///
    /// The decomposition splits V_local(r) = short(r) − Z·e²·erf(r)/r, where
    /// `short(r) = V_local(r) + Z·e²·erf(r)/r` is bounded at r=0 because
    /// erf(r)/r → 2/√π as r→0. FT[−Z·e²·erf(r)/r] = −4π·Z·e²·exp(−G²/4)/G².
    ///
    /// This is mathematically equivalent to the bare-Coulomb decomposition
    /// V(r) = [V(r) + Z·e²/r] − Z·e²/r, but numerically more robust: the
    /// bare form has a 1/r divergence near r=0 that only the r² factor tames,
    /// producing very large (billions of eV) values at the first radial grid
    /// points for high-Z PPs. The erf form has a bounded short-range integrand
    /// everywhere. For the PPs currently in use (Si, Fe, C), both forms agree
    /// to machine precision once Simpson's rule is applied over the log mesh
    /// (see `tests/vloc_erf_consistency.rs`).
    ///
    /// Units: r in Å, g_norm in Å⁻¹, returns eV (potential in reciprocal
    /// space per unit cell). The erf Gaussian width is 1 Å (matching the G
    /// units); V_local(G) is invariant to the choice of Gaussian width
    /// because the decomposition is exact for any positive width.
    pub fn v_local_of_g(&self, g_norm: f64, omega: f64) -> f64 {
        let four_pi = 4.0 * std::f64::consts::PI;

        if g_norm < G_ZERO_THRESHOLD {
            // G = 0 case: ∫ r² [V_loc(r) + Z e²/r] dr
            // (Bare Coulomb subtraction here matches QE. The bracketed term is
            // short-ranged because V_loc(r) → −Z·e²/r as r→∞.)
            let integrand: Vec<f64> = self
                .r_grid
                .iter()
                .zip(self.v_local.iter())
                .map(|(&r, &v)| {
                    let v_short = v + self.z_valence * E2 / r.max(1e-20);
                    r * r * v_short
                })
                .collect();
            let integral = simpson_integrate(&integrand, &self.rab);
            four_pi / omega * integral
        } else {
            // G ≠ 0: erf-subtracted integrand (bounded everywhere).
            //   integrand[i] = r²·[V_loc(r) + Z·e²·erf(r)/r]·sin(Gr)/(Gr)
            //                ≡ [r·V_loc(r) + Z·e²·erf(r)]·sin(Gr)/G
            // At r=0: r·V_loc(r) = 0 and Z·e²·erf(0) = 0, so integrand(0) = 0.
            let two_over_sqrt_pi = 2.0 / std::f64::consts::PI.sqrt();
            let integrand: Vec<f64> = self
                .r_grid
                .iter()
                .zip(self.v_local.iter())
                .map(|(&r, &v)| {
                    // erf(r)/r with small-r limit 2/√π.
                    let erf_over_r = if r < 1e-20 {
                        two_over_sqrt_pi
                    } else {
                        puruspe::erf(r) / r
                    };
                    // Smooth short-range part: V_loc(r) + Z·e²·erf(r)/r
                    let v_short = v + self.z_valence * E2 * erf_over_r;
                    // sin(Gr)/(Gr) with series fallback for small Gr.
                    let gr = g_norm * r;
                    let sinc = if gr < 1e-10 { 1.0 - gr * gr / 6.0 } else { gr.sin() / gr };
                    r * r * v_short * sinc
                })
                .collect();
            let integral = simpson_integrate(&integrand, &self.rab);
            let g2 = g_norm * g_norm;
            four_pi / omega * integral - four_pi * self.z_valence * E2 * (-g2 / 4.0).exp() / (omega * g2)
        }
    }

    /// Check if atomic charge density is available (non-zero).
    pub fn has_rho_atom(&self) -> bool {
        self.rho_atom.iter().any(|&v| v.abs() > 1e-20)
    }
}
