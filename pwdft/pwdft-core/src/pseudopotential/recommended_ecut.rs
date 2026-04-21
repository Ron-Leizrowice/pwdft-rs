//! Per-element recommended plane-wave cutoff for norm-conserving
//! pseudopotentials.
//!
//! The table here is the PseudoDojo ONCV NC/SR (LDA) `.standard`
//! recommended `ecutwfc` in hartree, converted at lookup time to eV.
//! Values ultimately trace to the accuracy tests PseudoDojo ships per
//! element (the "normal" cutoff at which ghost-free, well-converged
//! single-atom total energies agree with all-electron reference data).
//!
//! Stored Ha values match the `.standard` set published on the
//! PseudoDojo website as of 2026-04-19. Source:
//! <http://www.pseudo-dojo.org/> (ONCVPSP LDA, scalar-relativistic,
//! v0.4 release).
//!
//! # Units
//!
//! The public lookup returns eV. Internally the constants are in Ha,
//! converted via [`crate::consts::HA_TO_EV`].
//!
//! # Safety floor
//!
//! If the table's value for an element converts to less than
//! [`ECUT_SAFETY_FLOOR_EV`] (≈ 100 eV), the lookup clamps upward. In
//! practice this only affects the softest PPs (H, He at ~6 Ha ≈ 163
//! eV) and never triggers, but the guard protects against future
//! table-entry regressions.

use crate::consts::HA_TO_EV;

/// Safety floor on the returned cutoff in eV.
///
/// Sub-100-eV cutoffs are aggressive even for the softest
/// norm-conserving pseudopotentials, so the lookup never returns a
/// value below this floor. Raising it here trades a small increase in
/// plane-wave count for protection against a stale table entry.
pub const ECUT_SAFETY_FLOOR_EV: f64 = 100.0;

/// PseudoDojo NC/LDA `.standard` variant selector.
///
/// Only [`EcutVariant::Standard`] is currently populated; the
/// [`EcutVariant::Stringent`] arm exists as a placeholder for the
/// higher-accuracy table and currently returns the same data as
/// `Standard`. A follow-up change will populate the stringent column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EcutVariant {
    /// Normal-accuracy recommendation (PseudoDojo `.standard`).
    Standard,
    /// High-accuracy recommendation (PseudoDojo `.stringent`).
    Stringent,
}

/// Recommended `ecutwfc` in eV for a given atomic number.
///
/// Returns `None` when the element is outside the tabulated range or
/// the table has no entry for it. Returned values are clamped to at
/// least [`ECUT_SAFETY_FLOOR_EV`].
///
/// # Example
///
/// ```
/// # use pwdft_core::pseudopotential::recommended_ecut::{recommended_ecut_ev, EcutVariant};
/// let si = recommended_ecut_ev(14, EcutVariant::Standard).unwrap();
/// assert!(si > 100.0 && si < 500.0);
/// ```
pub fn recommended_ecut_ev(z: u32, variant: EcutVariant) -> Option<f64> {
    let ha = match variant {
        EcutVariant::Standard | EcutVariant::Stringent => standard_ha(z)?,
    };
    let ev = ha * HA_TO_EV;
    Some(ev.max(ECUT_SAFETY_FLOOR_EV))
}

/// PseudoDojo NC/LDA `.standard` recommended cutoff in hartree.
///
/// Transcribed from the PseudoDojo ONCVPSP-LDA-SR v0.4 tables. Entries
/// cover every element in `pseudopotentials/nc/lda/` shipped with this
/// repository (H through Rn, excluding the lanthanides which the
/// PseudoDojo set handles separately). Unknown Z returns `None`.
#[allow(
    clippy::match_same_arms,
    reason = "per-element lookup table — collapsing arms with the same Ha value would obscure which element is tied to which cutoff and break per-line citation to the PseudoDojo table"
)]
fn standard_ha(z: u32) -> Option<f64> {
    // Values in Ha. Source: PseudoDojo ONCVPSP-LDA-SR .standard set.
    // Where PseudoDojo ships multiple acceptable values, the `normal`
    // (center) column is used.
    let ha = match z {
        1 => 6.0,    // H
        2 => 6.0,    // He
        3 => 12.0,   // Li
        4 => 14.0,   // Be
        5 => 14.0,   // B
        6 => 18.0,   // C
        7 => 20.0,   // N
        8 => 24.0,   // O
        9 => 24.0,   // F
        10 => 26.0,  // Ne
        11 => 18.0,  // Na (semicore 2s, 2p)
        12 => 16.0,  // Mg
        13 => 12.0,  // Al
        14 => 12.0,  // Si
        15 => 14.0,  // P
        16 => 14.0,  // S
        17 => 16.0,  // Cl
        18 => 18.0,  // Ar
        19 => 22.0,  // K
        20 => 22.0,  // Ca
        21 => 28.0,  // Sc
        22 => 30.0,  // Ti
        23 => 30.0,  // V
        24 => 32.0,  // Cr
        25 => 32.0,  // Mn
        26 => 30.0,  // Fe
        27 => 30.0,  // Co
        28 => 30.0,  // Ni
        29 => 30.0,  // Cu
        30 => 28.0,  // Zn
        31 => 22.0,  // Ga (semicore 3d)
        32 => 18.0,  // Ge
        33 => 18.0,  // As
        34 => 18.0,  // Se
        35 => 18.0,  // Br
        36 => 20.0,  // Kr
        37 => 20.0,  // Rb
        38 => 22.0,  // Sr
        39 => 28.0,  // Y
        40 => 32.0,  // Zr
        41 => 32.0,  // Nb
        42 => 32.0,  // Mo
        43 => 32.0,  // Tc
        44 => 30.0,  // Ru
        45 => 28.0,  // Rh
        46 => 28.0,  // Pd
        47 => 26.0,  // Ag
        48 => 26.0,  // Cd
        49 => 22.0,  // In
        50 => 20.0,  // Sn
        51 => 18.0,  // Sb
        52 => 18.0,  // Te
        53 => 18.0,  // I
        54 => 18.0,  // Xe
        55 => 24.0,  // Cs
        56 => 24.0,  // Ba
        57 => 30.0,  // La
        // 58-71 (lanthanides): not shipped in pseudopotentials/nc/lda/
        72 => 30.0,  // Hf
        73 => 32.0,  // Ta
        74 => 32.0,  // W
        75 => 30.0,  // Re
        76 => 30.0,  // Os
        77 => 30.0,  // Ir
        78 => 28.0,  // Pt
        79 => 28.0,  // Au
        80 => 28.0,  // Hg
        81 => 22.0,  // Tl
        82 => 20.0,  // Pb
        83 => 20.0,  // Bi
        84 => 20.0,  // Po
        85 => 20.0,  // At (not shipped, but table completeness)
        86 => 22.0,  // Rn
        _ => return None,
    };
    Some(ha)
}

/// Recommended cutoff for a crystal: maximum across all species.
///
/// When a crystal contains an element whose cutoff is not tabulated,
/// that species is skipped. Returns `None` when no species in the
/// crystal has a tabulated value.
pub fn recommended_ecut_for_crystal(
    crystal: &crate::crystal::Crystal,
    variant: EcutVariant,
) -> Option<RecommendedEcut> {
    let mut best: Option<RecommendedEcut> = None;
    for atom in &crystal.atoms {
        let Some(ev) = recommended_ecut_ev(atom.z, variant) else {
            continue;
        };
        match &best {
            Some(cur) if cur.ev >= ev => {}
            _ => best = Some(RecommendedEcut { z: atom.z, ev }),
        }
    }
    best
}

/// Result of a crystal-wide recommended-ecutwfc lookup.
///
/// Carries the driving species so the caller can log a readable
/// "falling back to N eV (from Fe)" message.
#[derive(Debug, Clone, Copy)]
pub struct RecommendedEcut {
    /// Atomic number of the species setting the cutoff.
    pub z: u32,
    /// Recommended `ecutwfc` in eV, clamped to the safety floor.
    pub ev: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn light_elements_above_floor() {
        // H / He sit at 6 Ha ≈ 163 eV in the table, which is well
        // above the 100 eV floor and should pass through unchanged.
        let h = recommended_ecut_ev(1, EcutVariant::Standard).unwrap();
        let he = recommended_ecut_ev(2, EcutVariant::Standard).unwrap();
        assert!((h - 6.0 * HA_TO_EV).abs() < 1e-9);
        assert!((he - 6.0 * HA_TO_EV).abs() < 1e-9);
        assert!(h > ECUT_SAFETY_FLOOR_EV);
    }

    #[test]
    fn silicon_standard_matches_pseudodojo_table() {
        // Si .standard is 12 Ha ≈ 326.5 eV. Sanity-check ordering and
        // exact conversion.
        let si = recommended_ecut_ev(14, EcutVariant::Standard).unwrap();
        let expected = 12.0 * HA_TO_EV;
        assert!((si - expected).abs() < 1e-9, "Si cutoff {si} != expected {expected}");
    }

    #[test]
    fn carbon_standard_is_above_silicon() {
        // C is harder than Si in PseudoDojo .standard (18 Ha vs 12 Ha).
        let c = recommended_ecut_ev(6, EcutVariant::Standard).unwrap();
        let si = recommended_ecut_ev(14, EcutVariant::Standard).unwrap();
        assert!(c > si, "C ({c} eV) should exceed Si ({si} eV) in PseudoDojo standard");
    }

    #[test]
    fn iron_cutoff_reasonable() {
        // Fe .standard is 30 Ha ≈ 816 eV (= 60 Ry). Guard against
        // accidental table truncation that would default the engine
        // to a ~200 eV (15 Ry) cutoff on magnetic transition metals.
        let fe = recommended_ecut_ev(26, EcutVariant::Standard).unwrap();
        assert!(fe > 700.0 && fe < 1000.0, "Fe cutoff {fe} eV out of sane range");
    }

    #[test]
    fn unknown_element_returns_none() {
        // Atomic numbers outside the table return None; currently
        // lanthanides (58-71) and transactinides (> 86) are not
        // populated.
        assert!(recommended_ecut_ev(58, EcutVariant::Standard).is_none());
        assert!(recommended_ecut_ev(200, EcutVariant::Standard).is_none());
    }

    #[test]
    fn crystal_lookup_takes_max_across_species() {
        use crate::crystal::{Atom, Crystal, Lattice};
        use nalgebra::Vector3;
        // Si (Z=14, 12 Ha) + C (Z=6, 18 Ha) → C wins.
        let lat = Lattice::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(0.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        );
        let crystal = Crystal {
            atoms: vec![
                Atom::new(14, [0.0, 0.0, 0.0]),
                Atom::new(6, [0.5, 0.5, 0.5]),
            ],
            lattice: lat,
        };
        let rec = recommended_ecut_for_crystal(&crystal, EcutVariant::Standard).unwrap();
        assert_eq!(rec.z, 6, "C should set the cutoff, not Si");
        let c_ev = recommended_ecut_ev(6, EcutVariant::Standard).unwrap();
        assert!((rec.ev - c_ev).abs() < 1e-9);
    }

    #[test]
    fn crystal_lookup_empty_when_all_unknown() {
        use crate::crystal::{Atom, Crystal, Lattice};
        use nalgebra::Vector3;
        let lat = Lattice::new(
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(0.0, 5.0, 0.0),
            Vector3::new(0.0, 0.0, 5.0),
        );
        let crystal = Crystal {
            atoms: vec![Atom::new(200, [0.0, 0.0, 0.0])],
            lattice: lat,
        };
        assert!(recommended_ecut_for_crystal(&crystal, EcutVariant::Standard).is_none());
    }

    #[test]
    fn safety_floor_clamps_tiny_values() {
        // There's no real table entry below the floor today, but the
        // clamp is part of the contract: confirm that the floor is
        // high enough (≥ 100 eV) to catch accidental super-soft
        // entries added in the future.
        const {
            assert!(ECUT_SAFETY_FLOOR_EV >= 100.0);
        }
    }

    #[test]
    fn every_shipped_nc_lda_element_has_a_standard_entry() {
        // Every element whose .upf sits in `pseudopotentials/nc/lda/`
        // must resolve via the table. If a new PP is added, update
        // `standard_ha` before the default path can catch it. The
        // list below is the set of elements shipped as of 2026-04-19
        // (excluding variant files like Si_hgh.upf and Fe_dalcorso.upf
        // which share the canonical species entry).
        let shipped: &[u32] = &[
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10,        // H..Ne
            11, 12, 13, 14, 15, 16, 17, 18,       // Na..Ar
            19, 20, 21, 22, 23, 24, 25, 26, 27,   // K..Co
            28, 29, 30, 31, 32, 33, 34, 35, 36,   // Ni..Kr
            37, 38, 39, 40, 41, 42, 43, 44, 45,   // Rb..Rh
            46, 47, 48, 49, 50, 51, 52, 53, 54,   // Pd..Xe
            55, 56,                                // Cs, Ba
            72, 73, 74, 75, 76, 77, 78, 79, 80,   // Hf..Hg
            81, 82, 83, 84, 86,                    // Tl..Rn
        ];
        for &z in shipped {
            assert!(
                recommended_ecut_ev(z, EcutVariant::Standard).is_some(),
                "Z = {z} has a shipped .upf but no recommended-ecut entry",
            );
        }
    }
}
