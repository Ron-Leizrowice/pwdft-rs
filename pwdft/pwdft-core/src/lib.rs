//! Plane-wave Density Functional Theory (DFT) solver.
//!
//! A self-consistent Kohn-Sham SCF loop on a plane-wave basis with
//! norm-conserving pseudopotentials. LDA (Perdew-Zunger) is implemented
//! today; GGA (PBE) and hybrid functionals are in progress. Apple Metal
//! GPU acceleration is available behind the `gpu` feature flag.
//!
//! Entry points:
//! - [`scf::run_scf`] — full self-consistent calculation.
//! - [`bandstructure::compute_band_structure`] — non-self-consistent eigenvalues along a k-path.
//!
//! Module groups:
//! - **Crystal & basis:** [`crystal`], [`basis`], [`kpoints`], [`atoms`].
//! - **Pseudopotentials:** [`pseudopotential`].
//! - **Potentials & Hamiltonian:** [`potential`], [`hamiltonian`].
//! - **SCF loop:** [`scf`] (driver, mixing, energy, density, smearing).
//! - **Numerics:** [`fft`], [`eigensolver`], [`numerics`], [`ewald`].
//! - **Symmetry:** [`symmetry`].
//!
//! Internal units are eV for energies, Å for lengths, e/Å³ for densities.
//! Conversion constants live in [`consts`]; Ry/Bohr appear only at the
//! pseudopotential parsing boundary.
#![warn(unused_results)]
#![cfg_attr(test, allow(clippy::unwrap_used, clippy::panic, unused_results))]

pub mod bandstructure;
pub mod basis;
pub mod calculation;
pub mod consts;
pub mod crystal;
pub mod eigensolver;
pub mod error;
pub mod ewald;
pub mod fft;
#[cfg(feature = "gpu")]
pub mod gpu;
pub mod hamiltonian;
pub mod kpoints;
pub mod numerics;
pub mod potential;
pub mod pseudopotential;
pub mod scf;
pub mod settings;
pub mod functionals;
pub mod symmetry;
