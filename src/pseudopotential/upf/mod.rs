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
pub fn parse(content: &str) -> Result<PseudopotentialData> {
    convert::parse_body(content)
}
