//! Text-level XML helpers for UPF parsing.
//!
//! UPF v2 files are shallow, whitespace-tolerant XML, so we parse them
//! with simple string searches rather than pulling in a full XML
//! dependency. All helpers here are `pub(super)` and used only by
//! [`super::convert`].

use crate::error::{PwdftError, Result};

/// Extract an XML attribute value: `attr_name="value"`.
pub(super) fn extract_attr<'a>(content: &'a str, attr_name: &str) -> Option<&'a str> {
    let pattern = format!("{attr_name}=\"");
    let start = content.find(&pattern)? + pattern.len();
    let end = content[start..].find('"')? + start;
    Some(&content[start..end])
}

/// Extract `angular_momentum` from a PP_BETA tag.
///
/// Returns `Err(PwdftError::InvalidInput)` when the value parses as a
/// negative integer (physically meaningless; would later propagate as an
/// out-of-bounds index into the spherical-harmonic tables inside
/// `NonlocalPotential::new`), and `Err(PwdftError::Parse)` when the tag
/// or attribute is missing or the integer is unparseable. Rejecting at
/// load time keeps the failure mode "bad UPF file" rather than an
/// internal assertion deep in the solver (UPFV, 2026-04-18).
///
/// Zero (s-channel) and positive integers parse unchanged.
pub(super) fn extract_beta_angular_momentum(content: &str, tag: &str) -> Result<i32> {
    let tag_start = content
        .find(&format!("<{tag}"))
        .ok_or_else(|| PwdftError::Parse(format!("missing tag <{tag}>")))?;
    let tag_end = content[tag_start..]
        .find('>')
        .ok_or_else(|| PwdftError::Parse(format!("malformed tag <{tag}>")))?
        + tag_start;
    let tag_content = &content[tag_start..tag_end];
    let am_str = extract_attr(tag_content, "angular_momentum").ok_or_else(|| {
        PwdftError::Parse(format!("missing angular_momentum attribute in {tag}"))
    })?;
    let l: i32 = am_str
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("angular_momentum parse error in {tag}: {e}")))?;
    if l < 0 {
        return Err(PwdftError::InvalidInput(format!(
            "{tag}: angular_momentum must be non-negative (got {l})"
        )));
    }
    Ok(l)
}

/// Extract a block of floating-point data between `<TAG ...>` and `</TAG>`.
pub(super) fn extract_data_block(content: &str, tag: &str, expected_size: usize) -> Result<Vec<f64>> {
    let open_tag = format!("<{tag}");
    let close_tag = format!("</{tag}>");

    let tag_pos = content
        .find(&open_tag)
        .ok_or_else(|| PwdftError::Parse(format!("missing tag <{tag}>")))?;

    // Find end of opening tag
    let data_start = content[tag_pos..]
        .find('>')
        .ok_or_else(|| PwdftError::Parse(format!("malformed tag <{tag}>")))?
        + tag_pos
        + 1;

    let data_end = content[data_start..]
        .find(&close_tag)
        .ok_or_else(|| PwdftError::Parse(format!("missing closing tag </{tag}>")))?
        + data_start;

    let data_str = &content[data_start..data_end];
    let values: Vec<f64> = data_str
        .split_whitespace()
        .map(|s| {
            s.parse::<f64>()
                .map_err(|e| PwdftError::Parse(format!("float parse error in {tag}: {e} ({s})")))
        })
        .collect::<Result<Vec<f64>>>()?;

    if values.len() != expected_size {
        return Err(PwdftError::Parse(format!(
            "{tag}: expected {expected_size} values, got {}",
            values.len()
        )));
    }

    Ok(values)
}
