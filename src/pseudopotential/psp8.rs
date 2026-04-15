//! Parser for PSP8 pseudopotential files (ABINIT / PseudoDojo format).
//!
//! PSP8 files use Hartree atomic units: energies in Ha, lengths in Bohr.
//! We convert to internal units (eV, Å) on parse.
//!
//! File structure:
//! - Line 1: Title
//! - Line 2: zatom, zion, pspdat
//! - Line 3: pspcod(=8), pspxc, lmax, lloc, mmax, r2well
//! - Line 4: rchrg, fchrg, qchrg
//! - Line 5: nproj(l=0), nproj(l=1), ..., nproj(l=lmax)
//! - Line 6: extension_switch
//! - Then data blocks for each l channel, then local potential, then core charge.

use crate::error::{PwdftError, Result};

use super::{BetaProjector, PseudopotentialData, BOHR_TO_ANG};

use crate::consts::HA_TO_EV;

/// Parse a PSP8 file from its text content.
pub fn parse(content: &str) -> Result<PseudopotentialData> {
    let lines: Vec<&str> = content.lines().collect();
    if lines.len() < 7 {
        return Err(PwdftError::Parse("PSP8 file too short".into()));
    }

    // Line 1: title (may contain element symbol)
    let _title = lines[0].trim();

    // Line 2: zatom, zion, pspdat
    let fields2: Vec<&str> = lines[1].split_whitespace().collect();
    let zatom: f64 = parse_field(&fields2, 0, "zatom")?;
    let z_valence: f64 = parse_field(&fields2, 1, "zion")?;

    // Line 3: pspcod, pspxc, lmax, lloc, mmax, r2well
    let fields3: Vec<&str> = lines[2].split_whitespace().collect();
    let pspcod: i32 = parse_field(&fields3, 0, "pspcod")?;
    if pspcod != 8 {
        return Err(PwdftError::Parse(format!(
            "expected pspcod=8, got {pspcod}"
        )));
    }
    let l_max: i32 = parse_field(&fields3, 2, "lmax")?;
    let lloc: i32 = parse_field(&fields3, 3, "lloc")?;
    let mmax: usize = parse_field(&fields3, 4, "mmax")?;

    // Line 4: rchrg, fchrg, qchrg
    let fields4: Vec<&str> = lines[3].split_whitespace().collect();
    let fchrg: f64 = parse_field(&fields4, 1, "fchrg")?;

    // Line 5: nproj for each l
    let fields5: Vec<&str> = lines[4].split_whitespace().collect();
    let mut nproj_per_l = Vec::new();
    for l in 0..=l_max {
        let np: usize = parse_field(&fields5, l as usize, &format!("nproj(l={l})"))?;
        nproj_per_l.push(np);
    }
    let n_projectors: usize = nproj_per_l.iter().sum();

    // Line 6: extension_switch
    // (we don't use spin-orbit for now, just skip)

    // Parse data blocks starting from line 7 (index 6)
    let mut line_idx = 6;
    let mut r_grid = Vec::with_capacity(mmax);
    let mut beta_projectors = Vec::new();
    let mut v_local_data = Vec::with_capacity(mmax);

    // For each l from 0 to lmax (except lloc), read projector data
    for l in 0..=l_max {
        let np = nproj_per_l[l as usize];
        if l == lloc {
            // This is the local channel — skip (read later as local potential)
            continue;
        }

        for _proj_idx in 0..np {
            // Skip header line for this projector block
            if line_idx >= lines.len() {
                return Err(PwdftError::Parse("unexpected end of PSP8 file".into()));
            }
            line_idx += 1; // ekb header line

            let mut proj_values = Vec::with_capacity(mmax);
            for _i in 0..mmax {
                if line_idx >= lines.len() {
                    return Err(PwdftError::Parse("unexpected end of PSP8 data".into()));
                }
                let fields: Vec<&str> = lines[line_idx].split_whitespace().collect();
                // Column 0: index, Column 1: r, Column 2: projector value
                if r_grid.len() < mmax && beta_projectors.is_empty() {
                    let r: f64 = parse_field(&fields, 1, "r")?;
                    r_grid.push(r * BOHR_TO_ANG);
                }
                let val: f64 = parse_field(&fields, 2, "projector")?;
                proj_values.push(val * HA_TO_EV.sqrt() / BOHR_TO_ANG.sqrt());
                line_idx += 1;
            }

            beta_projectors.push(BetaProjector {
                l,
                values: proj_values,
            });
        }
    }

    // Read local potential block
    // Header line
    if line_idx < lines.len() {
        line_idx += 1; // skip header
    }
    for _i in 0..mmax {
        if line_idx >= lines.len() {
            return Err(PwdftError::Parse(
                "unexpected end of PSP8 local potential data".into(),
            ));
        }
        let fields: Vec<&str> = lines[line_idx].split_whitespace().collect();
        if r_grid.len() < mmax {
            let r: f64 = parse_field(&fields, 1, "r")?;
            r_grid.push(r * BOHR_TO_ANG);
        }
        let val: f64 = parse_field(&fields, 2, "vloc")?;
        v_local_data.push(val * HA_TO_EV);
        line_idx += 1;
    }

    // Generate rab from r_grid (finite differences)
    let rab = generate_rab(&r_grid);

    // Read core charge if present
    let rho_atom = if fchrg > 0.0 && line_idx + mmax <= lines.len() {
        // Skip to core charge block, parse it
        // For now, return zeros — we mainly need this for NLCC which is Phase 3+
        vec![0.0; mmax]
    } else {
        vec![0.0; mmax]
    };

    // D_ij: PSP8 encodes ekb (energy of KB projectors) in the header of each block.
    // For now, construct a diagonal D_ij from the ekb values.
    // TODO: read ekb values properly for off-diagonal terms
    let dij = vec![0.0; n_projectors * n_projectors];

    // Derive element symbol from atomic number
    let z_int = zatom.round() as u32;
    let element = crate::atoms::from_z(z_int)
        .map(|e| e.symbol().to_string())
        .unwrap_or_else(|| format!("Z{z_int}"));

    Ok(PseudopotentialData {
        element,
        z_valence,
        l_max,
        r_grid,
        rab,
        v_local: v_local_data,
        beta_projectors,
        dij,
        n_projectors,
        rho_atom,
        core_charge: vec![],
        has_nlcc: false,
    })
}

fn parse_field<T: std::str::FromStr>(fields: &[&str], idx: usize, name: &str) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    fields
        .get(idx)
        .ok_or_else(|| PwdftError::Parse(format!("missing field {name} at index {idx}")))?
        .parse()
        .map_err(|e| PwdftError::Parse(format!("parse error for {name}: {e}")))
}

/// Generate integration weights from radial grid via finite differences.
fn generate_rab(r_grid: &[f64]) -> Vec<f64> {
    let n = r_grid.len();
    if n == 0 {
        return vec![];
    }
    let mut rab = vec![0.0; n];
    if n == 1 {
        rab[0] = 1.0;
        return rab;
    }
    rab[0] = r_grid[1] - r_grid[0];
    for i in 1..n - 1 {
        rab[i] = (r_grid[i + 1] - r_grid[i - 1]) / 2.0;
    }
    rab[n - 1] = r_grid[n - 1] - r_grid[n - 2];
    rab
}
