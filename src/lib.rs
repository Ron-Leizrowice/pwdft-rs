#![cfg_attr(
    test,
    allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        reason = "ERR2 § Phase 0: in-src test modules are allowed to panic; CLAU consolidated 32 per-module allows to this single crate-level cfg_attr"
    )
)]

pub mod atoms;
pub mod bandstructure;
pub mod basis;
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
pub mod symmetry;
