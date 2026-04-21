use std::{collections::HashMap, path::PathBuf};

use clap::Parser;
use elements_rs::Element;
use itertools::Itertools;
use log::info;
use pwdft_core::{
    bandstructure,
    basis::BasisSet,
    kpoints,
    pseudopotential::UpfPseudoPotential,
    scf,
    settings::{InputSettings, KPointSettings},
};

#[derive(Parser)]
#[command(name = "pwdft-core", about = "Plane-wave DFT solver")]
struct Cli {
    /// Path to input file (YAML or TOML).
    #[arg(short, long)]
    input: PathBuf,

    /// Output file for band structure TSV (default: stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() -> pwdft_core::error::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let settings = load_settings(&cli.input)?;
    let crystal = settings.to_crystal()?;

    info!(
        "Crystal: {} atoms, lattice volume = {:.3} ų",
        crystal.atoms.len(),
        crystal.lattice.volume()
    );

    let basis = BasisSet::new(&crystal.lattice, ecut);
    info!("Basis set: {} plane waves at ecut = {ecut} eV", basis.len());

    let n_bands_fallback = basis.len().min(20);

    match &settings.kpoints {
        KPointSettings::BandPath { npoints, .. } => {
            let path_points = settings.to_high_sym_path().ok_or_else(|| {
                pwdft_core::error::PwdftError::InvalidInput("band_path mode requires a band path definition".into())
            })?;
            let (kpts, distances) = kpoints::high_symmetry_path(&path_points, *npoints, &crystal.lattice);

            let n_bands = settings.scf.n_bands.unwrap_or(n_bands_fallback);
            info!("Band structure: {} k-points, {n_bands} bands", kpts.len());

            let bs = bandstructure::compute_band_structure(&basis, &kpts, &distances, n_bands)?;

            match &cli.output {
                Some(path) => {
                    let mut file = std::fs::File::create(path)?;
                    bs.write_tsv(&mut file)?;
                    info!("Band structure written to {}", path.display());
                },
                None => {
                    let mut stdout = std::io::stdout().lock();
                    bs.write_tsv(&mut stdout)?;
                },
            }
        },
        KPointSettings::MonkhorstPack { grid, shift } => {
            let full_kpts = kpoints::monkhorst_pack(grid[0], grid[1], grid[2], *shift, &crystal.lattice);

            // Always reduce k-points via `reduce_kpoints`. When the user
            // disables symmetry, `settings.to_symmetry_info` returns the
            // identity-only group (no time reversal); in that case the
            // reduction is a no-op — every orbit has size 1 and the output
            // matches `monkhorst_pack` bit-for-bit (see
            // `reduce_kpoints_with_identity_only_preserves_full_grid`).
            //
            // A previous version short-circuited via `is_trivial()` when
            // `n_ops == 1`, but that also triggered for real P1 crystals
            // where `n_ops == 1` yet `has_time_reversal == true` — skipping
            // the k ↔ −k folding that physics demands. The uniform call
            // below avoids that regression and matches the treatment in
            // `symmetrize_density_g`, which already always-calls.
            let symmetry_info = settings.to_symmetry_info(&crystal);
            info!("Symmetry: {} space group operations", symmetry_info.n_ops);
            let kpts = pwdft_core::symmetry::kpoints::reduce_kpoints(
                &full_kpts,
                *grid,
                *shift,
                &symmetry_info,
                &crystal.lattice,
            );
            info!(
                "Monkhorst-Pack grid: {}×{}×{} shift={:?} = {} → {} IBZ k-points",
                grid[0],
                grid[1],
                grid[2],
                shift,
                full_kpts.len(),
                kpts.len()
            );

            // Load pseudopotentials
            let pseudopotentials: HashMap<Element, UpfPseudoPotential> = crystal
                .atoms
                .iter()
                .map(|a| a.symbol)
                .unique()
                .map(|sym| UpfPseudoPotential::load(sym).map(|pp| (sym, pp)))
                .collect::<Result<_, _>>()?;

            let params = settings.to_scf_params(n_bands_fallback);

            let result = scf::run_scf(&crystal, &basis, &kpts, &pseudopotentials, &params, &symmetry_info)?;

            info!("SCF converged in {} iterations", result.n_iterations);
            info!("Total energy: {:.6} eV", result.total_energy);
            info!("Fermi energy: {:.6} eV", result.fermi_energy);
            for (ik, evs) in result.eigenvalues.iter().enumerate() {
                if ik < 3 || ik == result.eigenvalues.len() - 1 {
                    info!(
                        "  k-point {ik}: bands = {:?}",
                        evs.iter().map(|e| format!("{e:.4}")).collect::<Vec<_>>()
                    );
                }
            }
        },
    }

    Ok(())
}

/// Load settings from a YAML input file.
fn load_settings(path: &std::path::Path) -> pwdft_core::error::Result<InputSettings> {
    InputSettings::from_yaml_file(path)
}
