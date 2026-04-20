//! Parser for UPF v2 pseudopotential files (Quantum ESPRESSO format).
//!
//! UPF files use Rydberg atomic units: energies in Ry, lengths in Bohr.
//! We convert to internal units (eV, Å) on parse.
//!
//! Layout:
//! - `xml` — text-level helpers that pull attribute values and numeric
//!   data blocks out of the UPF XML.
//! - `convert` — unit-conversion body that assembles a
//!   [`crate::pseudopotential::PseudopotentialData`] from the parsed text.
//!
//! Only [`parse`] is public outside this folder; the helpers are
//! `pub(super)` and must not leak.

mod convert;
mod xml;

use crate::error::Result;

use super::PseudopotentialData;

/// Parse a UPF v2 file from its text content.
///
/// # Errors
///
/// Returns [`crate::error::PwdftError::Parse`] when the UPF text cannot be
/// decoded into [`PseudopotentialData`]. Concrete
/// triggers (all surfaced as `PwdftError::Parse` with an explanatory
/// message):
/// - Missing or malformed required tags (`PP_HEADER`, `PP_MESH`, `PP_LOCAL`,
///   `PP_NONLOCAL`, `PP_BETA`, `PP_DIJ`, `PP_R`, `PP_RAB`).
/// - Mismatched array lengths (e.g. `PP_LOCAL` length differs from
///   `mesh_size`, or a `PP_BETA` projector has a different length than the
///   radial grid).
/// - Non-numeric or non-finite values inside a numeric block.
/// - Unsupported UPF version (only v2 text is handled).
pub fn parse(content: &str) -> Result<PseudopotentialData> {
    convert::parse_body(content)
}
