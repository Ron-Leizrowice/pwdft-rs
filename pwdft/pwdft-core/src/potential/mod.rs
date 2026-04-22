//! Potential terms of the Kohn-Sham Hamiltonian.
//!
//! The three pieces that, together with the kinetic operator, make up
//! the effective single-particle Hamiltonian:
//!
//! - [`local`] — ionic local potential V_local(r), assembled from the pseudopotential's V_local(G)
//!   via Bessel transform.
//! - [`nonlocal`] — Kleinman-Bylander separable non-local projectors, arbitrary angular momentum l
//!   via spherical-harmonic recurrence.
//! - [`xc`] — exchange-correlation functionals (LDA today; GGA, hybrid in progress) with dispatch
//!   over spin and functional choice.
//!
//! The Hartree term is assembled inline in the SCF driver (it depends
//! only on ρ(G) via a simple 1/|G|² kernel and does not need its own
//! module). The SCF pipeline is the caller that glues these together.

pub mod local;
pub mod nonlocal;
pub mod xc;
