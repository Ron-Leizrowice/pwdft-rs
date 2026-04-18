//! SCF progress reporting: per-iteration log lines and final-convergence
//! energy summary.
//!
//! Extracted from `scf::driver` and `scf::driver_spin` so that adding a
//! new metric touches one file instead of two. Non-spin and spin drivers
//! feed the same [`IterationReport`] struct, with spin-specific fields
//! (`delta_up`, `delta_down`, `magnetization`) carried in the optional
//! [`SpinIterationFields`].

use indicatif::ProgressBar;
use log::info;

use super::energy::EnergyComponents;

/// Extra fields emitted per-iteration by the spin-polarized driver.
pub(super) struct SpinIterationFields {
    pub delta_up: f64,
    pub delta_down: f64,
    pub magnetization: f64,
}

/// Per-iteration SCF progress record.
///
/// All energies in eV. `iter` is 0-indexed (display code adds 1).
pub(super) struct IterationReport {
    pub iter: usize,
    pub e_total: f64,
    pub e_harris: f64,
    pub hf_diff: f64,
    /// `None` on the first iteration (no previous energy).
    pub de: Option<f64>,
    pub delta: f64,
    /// Spin-specific fields; `None` for nspin=1.
    pub spin: Option<SpinIterationFields>,
}

/// Emit one progress-bar update and one `info!` line for an SCF iteration.
pub(super) fn log_iteration(pb: &ProgressBar, r: &IterationReport) {
    pb.set_position((r.iter + 1) as u64);
    match &r.spin {
        None => {
            pb.set_message(format!(
                "E={:.4} eV  Δρ={:.1e}",
                r.e_total, r.delta
            ));
            info!(
                "SCF iter {:>3}: E_KS={:.6} eV  E_HF={:.6} eV  |HF-KS|={:.2e}  dE={:>10}  Δρ={:.2e}",
                r.iter + 1,
                r.e_total,
                r.e_harris,
                r.hf_diff,
                r.de.map_or("N/A".to_string(), |de| format!("{de:.2e}")),
                r.delta
            );
        }
        Some(sp) => {
            pb.set_message(format!(
                "E={:.4} eV  Δρ={:.1e}  M={:.2} μB",
                r.e_total, r.delta, sp.magnetization
            ));
            info!(
                "SCF iter {:>3}: E_KS={:.6} eV  E_HF={:.6} eV  |HF-KS|={:.2e}  dE={:>10}  Δρ={:.2e} (↑{:.2e} ↓{:.2e})  M={:.3} μB",
                r.iter + 1,
                r.e_total,
                r.e_harris,
                r.hf_diff,
                r.de.map_or("N/A".to_string(), |de| format!("{de:.2e}")),
                r.delta,
                sp.delta_up,
                sp.delta_down,
                sp.magnetization
            );
        }
    }
}

/// Emit the converged-SCF energy summary: total energy, Harris-Foulkes,
/// free energy, σ→0 energy.
///
/// The non-spin driver additionally calls [`log_entropy`] after this; the
/// spin driver historically did not log entropy (see
/// `scf::driver::run_scf_unpolarized` vs `scf::driver_spin::run_scf_spin`).
pub(super) fn log_convergence_summary(
    e_total: f64,
    e_harris: f64,
    hf_diff: f64,
    free_energy: f64,
    energy_sigma0: f64,
) {
    info!("Energy (E_KS):   {e_total:.6} eV");
    info!("Harris-Foulkes:  {e_harris:.6} eV  (|HF-KS|={hf_diff:.2e})");
    info!("Free energy (F): {free_energy:.6} eV");
    info!("E sigma→0 (E₀):  {energy_sigma0:.6} eV");
}

/// Log entropy contribution if nonzero. Called by the non-spin driver only
/// (historical behavior — see `log_convergence_summary` doc).
pub(super) fn log_entropy(ts: f64, n_atoms: usize) {
    if ts.abs() > 1e-8 {
        info!(
            "Entropy (-TS):   {:.6} eV ({:.3} meV/atom)",
            -ts,
            -ts * 1000.0 / n_atoms as f64
        );
    }
}

/// Emit the per-component energy breakdown (VGC5 diagnostic) at SCF
/// convergence.
pub(super) fn log_components(c: &EnergyComponents, e_total: f64, n_electrons: f64) {
    let e_sum = c.e_kinetic
        + c.e_local
        + c.e_local_g0_shift
        + c.e_nonlocal
        + c.e_hartree
        + c.e_xc
        + c.e_ewald;
    info!("--- Per-component energies (eV) ---");
    info!("  E_band       = {:.6}", c.e_band);
    info!("  E_kinetic    = {:.6}", c.e_kinetic);
    info!("  E_local      = {:.6}", c.e_local);
    info!(
        "  E_local(G=0) = {:.6}  (= V_loc(G=0)·N_el, N_el={:.3})",
        c.e_local_g0_shift, n_electrons
    );
    info!("  E_nonlocal   = {:.6}", c.e_nonlocal);
    info!("  E_hartree    = {:.6}", c.e_hartree);
    info!("  E_xc         = {:.6}", c.e_xc);
    info!("  E_ewald      = {:.6}", c.e_ewald);
    info!(
        "  E_sum(comp)  = {e_sum:.6}   (vs E_KS {e_total:.6}, Δ={:.2e})",
        e_sum - e_total
    );
}
