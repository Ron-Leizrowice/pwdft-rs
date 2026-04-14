use std::path::PathBuf;

use clap::Parser;
use log::info;

use pwdft_rs::{
    bandstructure,
    basis::BasisSet,
    input::{InputFile, KPointsConfig},
    kpoints,
    scf,
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
    info!(
        "Basis set: {} plane waves at ecut = {} eV",
        basis.len(),
        config.system.ecut
    );

    let n_bands = config.system.n_bands.unwrap_or_else(|| basis.len().min(20));

    match &config.kpoints {
        KPointsConfig::BandPath { npoints, .. } => {
            let path_points = config.to_high_sym_path().unwrap();
            let (kpts, distances) =
                kpoints::high_symmetry_path(&path_points, *npoints, &crystal.lattice);

            info!("Band structure: {} k-points, {} bands", kpts.len(), n_bands);

            // If SCF config is present, run SCF first, then compute bands with converged potential
            // For now, compute free-electron bands
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
            let full_kpts = kpoints::monkhorst_pack(grid[0], grid[1], grid[2], &crystal.lattice);

            // Detect symmetry and reduce k-points
            let symmetry = pwdft_rs::symmetry::SymmetryInfo::from_crystal(&crystal, 1e-5);
            info!("Symmetry: {} space group operations", symmetry.n_ops);

            let kpts = pwdft_rs::symmetry::kpoints::reduce_kpoints(
                &full_kpts,
                *grid,
                &symmetry,
                &crystal.lattice,
            );
            info!(
                "Monkhorst-Pack grid: {}×{}×{} = {} → {} IBZ k-points",
                grid[0], grid[1], grid[2], full_kpts.len(), kpts.len()
            );

            // Load pseudopotentials
            let scf_config = config.scf.as_ref().expect(
                "SCF calculation requires [scf] section with pseudopotential paths in input file",
            );

            let input_dir = cli.input.parent().unwrap_or(std::path::Path::new("."));
            let mut pp_data = Vec::new();
            for atom_input in &config.system.atoms {
                let sym = &atom_input.symbol;
                if pp_data.iter().any(|pp: &pwdft_rs::pseudopotential::PseudopotentialData| pp.element == *sym) {
                    continue;
                }
                let pp_path = scf_config
                    .pseudopotentials
                    .get(sym)
                    .unwrap_or_else(|| panic!("no pseudopotential path for element {sym}"));
                let pp_path = input_dir.join(pp_path);
                info!("Loading pseudopotential for {sym}: {}", pp_path.display());
                let pp = pwdft_rs::pseudopotential::load(&pp_path)?;
                pp_data.push(pp);
            }
            let pp_refs: Vec<&pwdft_rs::pseudopotential::PseudopotentialData> =
                pp_data.iter().collect();

            let params = scf::ScfParams {
                n_bands,
                max_iter: scf_config.max_iter,
                conv_threshold: scf_config.conv_threshold,
                mixing_beta: scf_config.mixing_beta,
                mixing_ndim: scf_config.mixing_ndim,
                smearing_sigma: scf_config.smearing_sigma,
                ecutrho_ratio: scf_config.ecutrho_ratio,
                fft_grid: scf_config.fft_grid,
                mixing_mode: scf::mixing::MixingMode::Kerker { q_tf: None },
            };

            let result = scf::run_scf(&crystal, &basis, &kpts, &pp_refs, &params, Some(&symmetry))?;

            eprintln!("SCF converged in {} iterations", result.n_iterations);
            eprintln!("Total energy: {:.6} eV", result.total_energy);
            eprintln!("Fermi energy: {:.6} eV", result.fermi_energy);
            for (ik, evs) in result.eigenvalues.iter().enumerate() {
                if ik < 3 || ik == result.eigenvalues.len() - 1 {
                    eprintln!(
                        "  k-point {}: bands = {:?}",
                        ik,
                        evs.iter().map(|e| format!("{e:.4}")).collect::<Vec<_>>()
                    );
                }
            }
        }
    }

    Ok(())
}
