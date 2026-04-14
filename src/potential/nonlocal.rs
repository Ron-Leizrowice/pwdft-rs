//! Non-local pseudopotential in the Kleinman-Bylander separable form.
//!
//! V_NL = Σ_{atom} Σ_{i,j} |β_i⟩ D_{ij} ⟨β_j| × S_atom
//!
//! In reciprocal space at k-point k:
//! V_NL(k)_{G,G'} = Σ_{atom} S_atom(G-G') × Σ_{i,j} β_i(k+G) D_{ij} β_j*(k+G')
//!
//! where β_i(k+G) is the Fourier transform of the real-space projector:
//! β_{l}(q) = 4π i^l ∫ r·β_l(r) j_l(qr) r dr × Y_lm(q̂)
//!
//! (UPF stores r·β(r), so the integral is ∫ [r·β(r)] j_l(qr) r dr)

use nalgebra::Vector3;
use num_complex::Complex64;
use std::f64::consts::PI;

use crate::{
    basis::BasisSet,
    crystal::Crystal,
    pseudopotential::PseudopotentialData,
};

/// Precomputed non-local projector form factors β_i(q) for all |q| values needed.
pub struct NonlocalPotential {
    /// For each atom type: form factors β_i(|k+G|) for each projector.
    /// Indexed as [atom_type][projector_index][g_index].
    form_factors: Vec<Vec<Vec<f64>>>,
    /// D_ij matrices for each atom type (n_proj × n_proj, row-major).
    dij: Vec<Vec<f64>>,
    /// Number of projectors per atom type.
    n_proj: Vec<usize>,
    /// Angular momentum of each projector per atom type.
    proj_l: Vec<Vec<i32>>,
}

impl NonlocalPotential {
    /// Precompute projector form factors for a given k-point.
    ///
    /// For each projector i with angular momentum l:
    /// F_i(|k+G|) = 4π ∫ [r·β_i(r)] j_l(|k+G|·r) r dr
    ///
    /// The full projector in G-space is:
    /// β_i(k+G) = F_i(|k+G|) × i^l × Y_lm(k̂+G)
    ///
    /// But for the KB matrix element, we sum over m:
    /// Σ_m Y_lm(q̂) Y_lm*(q̂') = (2l+1)/(4π) P_l(cos θ)
    /// where θ is the angle between q and q'.
    ///
    /// So: V_NL_{G,G'} = Σ_atom S(G-G') Σ_{i,j with same l}
    ///     F_i(|k+G|) D_{ij} F_j(|k+G'|) × (2l+1)/(4π) P_l(cos θ) × (-1)^l
    ///
    /// (The i^l × (i*)^l = (-1)^l factor)
    pub fn new(
        crystal: &Crystal,
        basis: &BasisSet,
        k: &Vector3<f64>,
        pseudopotentials: &[&PseudopotentialData],
    ) -> Self {
        // Identify unique atom types
        let mut atom_types: Vec<u32> = crystal.atoms.iter().map(|a| a.z).collect();
        atom_types.sort();
        atom_types.dedup();

        let n_pw = basis.len();
        let mut form_factors = Vec::new();
        let mut dij_all = Vec::new();
        let mut n_proj_all = Vec::new();
        let mut proj_l_all = Vec::new();

        for &z in &atom_types {
            let pp = crate::pseudopotential::find_for_atom(z, pseudopotentials);

            let mut type_ff = Vec::new();
            let mut type_l = Vec::new();

            for proj in &pp.beta_projectors {
                let l = proj.l;
                type_l.push(l);

                // Compute F(|k+G|) for each G-vector
                let mut ff = Vec::with_capacity(n_pw);
                for g in basis.g_vectors() {
                    let q = k + g;
                    let q_norm = q.norm();
                    let f = bessel_transform_projector(&pp.r_grid, &pp.rab, &proj.values, l, q_norm);
                    ff.push(f);
                }
                type_ff.push(ff);
            }

            form_factors.push(type_ff);
            dij_all.push(pp.dij.clone());
            n_proj_all.push(pp.n_projectors);
            proj_l_all.push(type_l);
        }

        Self {
            form_factors,
            dij: dij_all,
            n_proj: n_proj_all,
            proj_l: proj_l_all,
        }
    }

    /// Add V_NL to the Hamiltonian matrix at a given k-point.
    ///
    /// V_NL(G,G') = (1/Ω) Σ_atom S(G-G') × Σ_{i,j} F_i(|k+G|) D_{ij} F_j(|k+G'|)
    ///              × (2l+1)/(4π) P_l(cos θ)
    ///
    /// Optimized: structure factors factored as exp(-iG·τ) × exp(iG'·τ),
    /// q-vectors and norms precomputed, projector sums lifted out of atom loop.
    pub fn add_to_hamiltonian(
        &self,
        h: &mut nalgebra::DMatrix<Complex64>,
        crystal: &Crystal,
        basis: &BasisSet,
        k: &Vector3<f64>,
    ) {
        let n_pw = basis.len();
        let omega = crystal.lattice.volume().abs();
        let inv_omega = 1.0 / omega;

        // Precompute q-vectors and norms (once, not per pair)
        let g_vecs = basis.g_vectors();
        let q_vecs: Vec<Vector3<f64>> = g_vecs.iter().map(|g| k + g).collect();
        let q_norms: Vec<f64> = q_vecs.iter().map(|q| q.norm()).collect();

        // Identify unique atom types
        let mut atom_types: Vec<u32> = crystal.atoms.iter().map(|a| a.z).collect();
        atom_types.sort();
        atom_types.dedup();

        for (itype, &z) in atom_types.iter().enumerate() {
            let n_proj = self.n_proj[itype];
            let dij = &self.dij[itype];
            let proj_l = &self.proj_l[itype];

            // Precompute per-atom structure factor phases: exp(-iG·τ) for each G
            let atoms_of_type: Vec<_> = crystal.atoms.iter().filter(|a| a.z == z).collect();
            let atom_phases: Vec<Vec<Complex64>> = atoms_of_type
                .iter()
                .map(|atom| {
                    let tau = atom.cart_position(&crystal.lattice);
                    g_vecs
                        .iter()
                        .map(|g| {
                            let phase = -g.dot(&tau);
                            Complex64::new(phase.cos(), phase.sin())
                        })
                        .collect()
                })
                .collect();

            for ig in 0..n_pw {
                let q_i = &q_vecs[ig];
                let q_i_norm = q_norms[ig];

                for jg in 0..n_pw {
                    let q_j = &q_vecs[jg];
                    let q_j_norm = q_norms[jg];

                    // Projector sum (same for all atoms of this type)
                    let mut vnl = 0.0;
                    for i in 0..n_proj {
                        for j in 0..n_proj {
                            if proj_l[i] != proj_l[j] {
                                continue;
                            }
                            let l = proj_l[i];

                            let fi = self.form_factors[itype][i][ig];
                            let fj = self.form_factors[itype][j][jg];
                            let d = dij[i * n_proj + j];

                            let cos_theta = if q_i_norm > 1e-12 && q_j_norm > 1e-12 {
                                q_i.dot(q_j) / (q_i_norm * q_j_norm)
                            } else {
                                1.0
                            };
                            let angular =
                                (2 * l + 1) as f64 / (4.0 * PI) * legendre_p(l, cos_theta);

                            vnl += fi * d * fj * angular;
                        }
                    }

                    if vnl.abs() < 1e-20 {
                        continue;
                    }

                    // Sum structure factors over atoms: Σ_atom exp(-iG_i·τ) × exp(iG_j·τ)
                    let mut sf_sum = Complex64::new(0.0, 0.0);
                    for phases in &atom_phases {
                        // S(G_i - G_j) = exp(-iG_i·τ) × conj(exp(-iG_j·τ))
                        sf_sum += phases[ig] * phases[jg].conj();
                    }

                    h[(ig, jg)] += sf_sum * (vnl * inv_omega);
                }
            }
        }
    }
}

/// Spherical Bessel transform of a projector:
/// F(q) = 4π ∫₀^∞ [r·β(r)] j_l(qr) r dr
///
/// where `r_beta` stores r·β(r) (the UPF convention).
fn bessel_transform_projector(
    r_grid: &[f64],
    rab: &[f64],
    r_beta: &[f64],
    l: i32,
    q: f64,
) -> f64 {
    let mut integral = 0.0;

    for i in 0..r_grid.len() {
        let r = r_grid[i];
        let dr = rab[i];
        let rb = r_beta[i]; // r·β(r)

        let qr = q * r;
        let jl = spherical_bessel_j(l, qr);

        // Integrand: [r·β(r)] × j_l(qr) × r × dr
        integral += rb * jl * r * dr;
    }

    4.0 * PI * integral
}

/// Spherical Bessel function j_l(x) for l = 0, 1, 2, 3.
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
        _ => panic!("spherical_bessel_j: l={l} not implemented (max l=3)"),
    }
}

/// Legendre polynomial P_l(x) for l = 0, 1, 2, 3.
fn legendre_p(l: i32, x: f64) -> f64 {
    match l {
        0 => 1.0,
        1 => x,
        2 => 0.5 * (3.0 * x * x - 1.0),
        3 => 0.5 * (5.0 * x * x * x - 3.0 * x),
        _ => panic!("legendre_p: l={l} not implemented (max l=3)"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[test]
    fn test_spherical_bessel_j0() {
        assert!(relative_eq!(spherical_bessel_j(0, 0.0), 1.0, epsilon = 1e-10));
        // j_0(π) = sin(π)/π = 0
        assert!(spherical_bessel_j(0, PI).abs() < 1e-10);
        // j_0(1) = sin(1)/1 ≈ 0.8415
        assert!(relative_eq!(
            spherical_bessel_j(0, 1.0),
            1.0_f64.sin(),
            epsilon = 1e-10
        ));
    }

    #[test]
    fn test_spherical_bessel_j1() {
        assert!(relative_eq!(spherical_bessel_j(1, 0.0), 0.0, epsilon = 1e-10));
        // j_1(x) = sin(x)/x² - cos(x)/x
        let x: f64 = 2.0;
        let expected = x.sin() / (x * x) - x.cos() / x;
        assert!(relative_eq!(
            spherical_bessel_j(1, x),
            expected,
            epsilon = 1e-10
        ));
    }

    #[test]
    fn test_legendre_orthogonality() {
        // ∫₋₁¹ P_l(x) P_m(x) dx = 2/(2l+1) δ_{lm}
        // Approximate with Gauss quadrature (Simpson's rule on [-1,1])
        let n = 1000;
        for l in 0..=3 {
            for m in 0..=3 {
                let mut integral = 0.0;
                for i in 0..=n {
                    let x = -1.0 + 2.0 * i as f64 / n as f64;
                    let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
                    integral += w * legendre_p(l, x) * legendre_p(m, x);
                }
                integral *= 2.0 / (3.0 * n as f64);
                let expected = if l == m {
                    2.0 / (2 * l + 1) as f64
                } else {
                    0.0
                };
                assert!(
                    (integral - expected).abs() < 0.01,
                    "P_{l} · P_{m} integral: {integral:.6}, expected {expected:.6}"
                );
            }
        }
    }

    #[test]
    fn test_bessel_transform_delta() {
        // For a delta-like projector at r=0, the Bessel transform is constant
        // Not easily testable, but ensure it returns finite values
        let r = vec![0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0];
        let rab = vec![0.01, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0];
        let rbeta: Vec<f64> = r.iter().map(|&ri: &f64| (-ri * ri).exp()).collect();

        let f0 = bessel_transform_projector(&r, &rab, &rbeta, 0, 0.0);
        let f1 = bessel_transform_projector(&r, &rab, &rbeta, 0, 1.0);
        let f5 = bessel_transform_projector(&r, &rab, &rbeta, 0, 5.0);

        assert!(f0.is_finite());
        assert!(f1.is_finite());
        assert!(f5.is_finite());
        // Should decay with q
        assert!(f5.abs() < f0.abs());
    }
}
