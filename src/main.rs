use std::path::PathBuf;

use clap::Parser;
use log::info;

use pwdft_rs::{
    bandstructure,
    basis::BasisSet,
    kpoints,
    scf,
    settings::{KPointSettings, Settings},
};

#[derive(Parser)]
#[command(name = "pwdft-rs", about = "Plane-wave DFT solver")]
struct Cli {
    /// Path to input file (YAML or TOML).
    #[arg(short, long)]
    input: PathBuf,

    /// Output file for band structure TSV (default: stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() -> pwdft_rs::error::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let settings = load_settings(&cli.input)?;
    let crystal = settings.to_crystal();

    info!(
        "Crystal: {} atoms, lattice volume = {:.3} ų",
        crystal.atoms.len(),
        crystal.lattice.volume()
    );

    let ecut = settings.ecutwfc();
    let basis = BasisSet::new(&crystal.lattice, ecut);
    info!("Basis set: {} plane waves at ecut = {ecut} eV", basis.len());

    let n_bands_fallback = basis.len().min(20);

    match &settings.kpoints {
        KPointSettings::BandPath { npoints, .. } => {
            let path_points = settings.to_high_sym_path().unwrap();
            let (kpts, distances) =
                kpoints::high_symmetry_path(&path_points, *npoints, &crystal.lattice);

            let n_bands = settings.scf.n_bands.unwrap_or(n_bands_fallback);
            info!("Band structure: {} k-points, {n_bands} bands", kpts.len());

            let bs =
                bandstructure::compute_band_structure(&basis, &kpts, &distances, n_bands, None);

            match &cli.output {
                Some(path) => {
                    let mut file = std::fs::File::create(path)?;
                    bs.write_tsv(&mut file)?;
                    info!("Band structure written to {}", path.display());
                }
                None => {
                    let mut stdout = std::io::stdout().lock();
                    bs.write_tsv(&mut stdout)?;
                }
            }
        }
        KPointSettings::MonkhorstPack { grid } => {
            let full_kpts =
                kpoints::monkhorst_pack(grid[0], grid[1], grid[2], &crystal.lattice);

            // Detect symmetry and reduce k-points
            let symmetry_info = settings.to_symmetry_info(&crystal);
            let kpts = if let Some(ref sym) = symmetry_info {
                info!("Symmetry: {} space group operations", sym.n_ops);
                let reduced = pwdft_rs::symmetry::kpoints::reduce_kpoints(
                    &full_kpts, *grid, sym, &crystal.lattice,
                );
                info!(
                    "Monkhorst-Pack grid: {}×{}×{} = {} → {} IBZ k-points",
                    grid[0], grid[1], grid[2], full_kpts.len(), reduced.len()
                );
                reduced
            } else {
                info!(
                    "Monkhorst-Pack grid: {}×{}×{} = {} k-points (no symmetry)",
                    grid[0], grid[1], grid[2], full_kpts.len()
                );
                full_kpts
            };

            // Load pseudopotentials
            let input_dir = cli.input.parent().unwrap_or(std::path::Path::new("."));
            let mut pp_data = Vec::new();
            for atom in &crystal.atoms {
                let sym = pwdft_rs::atoms::from_z(atom.z)
                    .map(|e| e.symbol().to_string())
                    .unwrap_or_default();
                if pp_data
                    .iter()
                    .any(|pp: &pwdft_rs::pseudopotential::PseudopotentialData| pp.element == sym)
                {
                    continue;
                }
                let pp_path = settings
                    .pseudopotential_path(&sym)
                    .unwrap_or_else(|| panic!("no pseudopotential path for element {sym}"));
                let pp_path = input_dir.join(pp_path);
                info!("Loading pseudopotential for {sym}: {}", pp_path.display());
                let pp = pwdft_rs::pseudopotential::load(&pp_path)?;
                pp_data.push(pp);
            }
            let pp_refs: Vec<&pwdft_rs::pseudopotential::PseudopotentialData> =
                pp_data.iter().collect();

            let params = settings.to_scf_params(n_bands_fallback);

            let result = scf::run_scf(
                &crystal,
                &basis,
                &kpts,
                &pp_refs,
                &params,
                symmetry_info.as_ref(),
            )?;

            eprintln!("SCF converged in {} iterations", result.n_iterations);
            eprintln!("Total energy: {:.6} eV", result.total_energy);
            eprintln!("Fermi energy: {:.6} eV", result.fermi_energy);
            for (ik, evs) in result.eigenvalues.iter().enumerate() {
                if ik < 3 || ik == result.eigenvalues.len() - 1 {
                    eprintln!(
                        "  k-point {ik}: bands = {:?}",
                        evs.iter().map(|e| format!("{e:.4}")).collect::<Vec<_>>()
                    );
                }
            }
        }
    }

    Ok(())
}

/// Load settings from a YAML input file.
fn load_settings(path: &std::path::Path) -> pwdft_rs::error::Result<Settings> {
    Settings::from_yaml_file(path)
}
