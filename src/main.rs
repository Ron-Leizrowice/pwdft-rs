use std::path::PathBuf;

use clap::Parser;
use log::info;

use pwdft_rs::{
    bandstructure,
    basis::BasisSet,
    input::{InputFile, KPointsConfig},
    kpoints,
};

#[derive(Parser)]
#[command(name = "pwdft-rs", about = "Plane-wave DFT solver")]
struct Cli {
    /// Path to TOML input file.
    #[arg(short, long)]
    input: PathBuf,

    /// Output file for band structure TSV (default: stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() -> pwdft_rs::error::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let config = InputFile::from_file(&cli.input)?;
    let crystal = config.to_crystal();

    info!(
        "Crystal: {} atoms, lattice volume = {:.3} ų",
        crystal.atoms.len(),
        crystal.lattice.volume()
    );

    let basis = BasisSet::new(&crystal.lattice, config.system.ecut);
    info!("Basis set: {} plane waves at ecut = {} eV", basis.len(), config.system.ecut);

    let n_bands = config.system.n_bands.unwrap_or_else(|| {
        // Default: enough bands to cover interesting physics
        (basis.len()).min(20)
    });

    match &config.kpoints {
        KPointsConfig::BandPath { npoints, .. } => {
            let path_points = config.to_high_sym_path().unwrap();
            let (kpts, distances) =
                kpoints::high_symmetry_path(&path_points, *npoints, &crystal.lattice);

            info!(
                "Band structure: {} k-points, {} bands",
                kpts.len(),
                n_bands
            );

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
        KPointsConfig::MonkhorstPack { grid } => {
            let kpts = kpoints::monkhorst_pack(grid[0], grid[1], grid[2], &crystal.lattice);
            info!(
                "Monkhorst-Pack grid: {}×{}×{} = {} k-points",
                grid[0],
                grid[1],
                grid[2],
                kpts.len()
            );
            eprintln!("SCF calculation not yet implemented. Use band_path kpoints for free-electron band structure.");
        }
    }

    Ok(())
}
