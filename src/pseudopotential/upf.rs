//! Parser for UPF v2 pseudopotential files (Quantum ESPRESSO format).
//!
//! UPF files use Rydberg atomic units: energies in Ry, lengths in Bohr.
//! We convert to internal units (eV, Å) on parse.

use crate::error::{PwdftError, Result};

use super::{BetaProjector, PseudopotentialData, BOHR_TO_ANG, RY_TO_EV};

/// Parse a UPF v2 file from its text content.
pub fn parse(content: &str) -> Result<PseudopotentialData> {
    let element = extract_attr(content, "element")
        .ok_or_else(|| PwdftError::Parse("missing element in PP_HEADER".into()))?
        .trim()
        .to_string();

    let z_valence: f64 = extract_attr(content, "z_valence")
        .ok_or_else(|| PwdftError::Parse("missing z_valence".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("z_valence parse error: {e}")))?;

    let l_max: i32 = extract_attr(content, "l_max")
        .ok_or_else(|| PwdftError::Parse("missing l_max".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("l_max parse error: {e}")))?;

    let n_proj: usize = extract_attr(content, "number_of_proj")
        .ok_or_else(|| PwdftError::Parse("missing number_of_proj".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("number_of_proj parse error: {e}")))?;

    let mesh_size: usize = extract_attr(content, "mesh_size")
        .ok_or_else(|| PwdftError::Parse("missing mesh_size".into()))?
        .trim()
        .parse()
        .map_err(|e| PwdftError::Parse(format!("mesh_size parse error: {e}")))?;

    // Parse radial grid (Bohr → Å)
    let r_grid_bohr = extract_data_block(content, "PP_R", mesh_size)?;
    let r_grid: Vec<f64> = r_grid_bohr.iter().map(|&r| r * BOHR_TO_ANG).collect();

    // Parse integration weights (Bohr → Å)
    let rab_bohr = extract_data_block(content, "PP_RAB", mesh_size)?;
    let rab: Vec<f64> = rab_bohr.iter().map(|&dr| dr * BOHR_TO_ANG).collect();

    // Parse local potential (Ry → eV)
    let v_local_ry = extract_data_block(content, "PP_LOCAL", mesh_size)?;
    let v_local: Vec<f64> = v_local_ry.iter().map(|&v| v * RY_TO_EV).collect();

    // Parse non-local projectors
    // UPF stores β(r) · r, in Ry^{1/2} units, on the radial grid.
    // We need β(r) in our units.
    let mut beta_projectors = Vec::with_capacity(n_proj);
    for i in 1..=n_proj {
        let tag = format!("PP_BETA.{i}");
        let l = extract_beta_angular_momentum(content, &tag)
            .ok_or_else(|| PwdftError::Parse(format!("missing angular_momentum in {tag}")))?;

        let beta_r_ry = extract_data_block(content, &tag, mesh_size)?;
        // UPF stores χ(r) = r · β(r) in Bohr^{-1/2} (no energy dimension).
        // β(r) is a wavefunction-like quantity in Bohr^{-3/2}, so χ = r·β is in Bohr^{-1/2}.
        // Energy enters only through D_ij (in Ry, converted to eV).
        //
        // Convert Bohr^{-1/2} → Å^{-1/2}: divide by √(BOHR_TO_ANG).
        let values: Vec<f64> = beta_r_ry
            .iter()
            .map(|&v| v / BOHR_TO_ANG.sqrt())
            .collect();

        beta_projectors.push(BetaProjector { l, values });
    }

    // Parse D_ij matrix (Ry → eV)
    let dij_ry = extract_data_block(content, "PP_DIJ", n_proj * n_proj)?;
    let dij: Vec<f64> = dij_ry.iter().map(|&d| d * RY_TO_EV).collect();

    // Parse atomic charge density (optional — may be all zeros for HGH)
    // UPF stores 4πr²ρ(r) in e/Bohr on the radial grid.
    // The Bessel transform ∫ [4πr²ρ(r)] j₀(Gr) dr needs consistent units:
    // r_grid is in Å, rab is in Å, G is in 1/Å, so 4πr²ρ(r) must be in e/Å.
    // Convert e/Bohr → e/Å by dividing by BOHR_TO_ANG.
    let rho_atom = if let Ok(rho_raw) = extract_data_block(content, "PP_RHOATOM", mesh_size) {
        rho_raw
            .iter()
            .map(|&rho| rho / BOHR_TO_ANG)
            .collect()
    } else {
        vec![0.0; mesh_size]
    };

    // Parse nonlinear core correction (NLCC) charge density.
    // PP_NLCC stores 4πr²ρ_core(r) in e/Bohr on the radial grid.
    // Same unit conversion as rho_atom: e/Bohr → e/Å.
    let has_nlcc = content.contains("core_correction=\"T\"")
        || content.contains("core_correction=\"t\"")
        || content.contains("nlcc=.true.");
    let core_charge = if has_nlcc {
        if let Ok(nlcc_raw) = extract_data_block(content, "PP_NLCC", mesh_size) {
            nlcc_raw
                .iter()
                .map(|&rho| rho / BOHR_TO_ANG)
                .collect()
        } else {
            vec![0.0; mesh_size]
        }
    } else {
        vec![]
    };

    Ok(PseudopotentialData {
        element,
        z_valence,
        l_max,
        r_grid,
        rab,
        v_local,
        beta_projectors,
        dij,
        n_projectors: n_proj,
        rho_atom,
        core_charge,
        has_nlcc,
    })
}

/// Extract an XML attribute value: `attr_name="value"`.
fn extract_attr<'a>(content: &'a str, attr_name: &str) -> Option<&'a str> {
    let pattern = format!("{attr_name}=\"");
    let start = content.find(&pattern)? + pattern.len();
    let end = content[start..].find('"')? + start;
    Some(&content[start..end])
}

/// Extract angular_momentum from a PP_BETA tag.
fn extract_beta_angular_momentum(content: &str, tag: &str) -> Option<i32> {
    // Find the tag opening
    let tag_start = content.find(&format!("<{tag}"))?;
    let tag_end = content[tag_start..].find('>')? + tag_start;
    let tag_content = &content[tag_start..tag_end];
    let am_str = extract_attr(tag_content, "angular_momentum")?;
    am_str.trim().parse().ok()
}

/// Extract a block of floating-point data between `<TAG ...>` and `</TAG>`.
fn extract_data_block(content: &str, tag: &str, expected_size: usize) -> Result<Vec<f64>> {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn si_content() -> String {
        std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/nc/lda/Si.upf"),
        )
        .unwrap()
    }

    #[test]
    fn test_parse_header() {
        let pp = parse(&si_content()).unwrap();
        assert_eq!(pp.element, "Si");
        assert!((pp.z_valence - 4.0).abs() < 1e-10);
        assert!(pp.l_max >= 1, "Si should have l_max >= 1");
        assert!(pp.n_projectors > 0, "Si should have projectors");
    }

    #[test]
    fn test_radial_grid_monotonic() {
        let pp = parse(&si_content()).unwrap();
        for i in 1..pp.r_grid.len() {
            assert!(
                pp.r_grid[i] > pp.r_grid[i - 1],
                "radial grid not monotonic at index {i}"
            );
        }
    }

    #[test]
    fn test_radial_grid_units() {
        let pp = parse(&si_content()).unwrap();
        // First grid point should be small (< 0.01 Å)
        assert!(pp.r_grid[0] < 0.01, "first r = {} Å too large", pp.r_grid[0]);
        // Last grid point should be reasonable (> 1 Å)
        assert!(
            pp.r_grid.last().unwrap() > &1.0,
            "last r = {} Å too small",
            pp.r_grid.last().unwrap()
        );
    }

    #[test]
    fn test_v_local_coulomb_tail() {
        let pp = parse(&si_content()).unwrap();
        // At large r, V_local should approach -Z_val e²/r
        // e² = 14.3997 eV·Å, Z_val = 4
        let e2 = 14.399645351950548;
        let n = pp.r_grid.len();
        let r_far = pp.r_grid[n - 10];
        let v_far = pp.v_local[n - 10];
        let v_coulomb = -pp.z_valence * e2 / r_far;
        let relative_err = ((v_far - v_coulomb) / v_coulomb).abs();
        assert!(
            relative_err < 0.01,
            "V_local at r={r_far:.2} Å: {v_far:.6} eV vs Coulomb {v_coulomb:.6} eV (err={relative_err:.4})"
        );
    }

    #[test]
    fn test_dij_matrix_size() {
        let pp = parse(&si_content()).unwrap();
        assert_eq!(pp.dij.len(), pp.n_projectors * pp.n_projectors);
    }
}
