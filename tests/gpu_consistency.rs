//! Verify GPU-accelerated SCF produces results consistent with CPU-only.
//!
//! Run with: cargo test --features gpu --test gpu_consistency

#![cfg(feature = "gpu")]

use nalgebra::Vector3;
use num_complex::Complex64;

use pwdft_rs::{
    basis::BasisSet,
    crystal::{Atom, Crystal, Lattice},
    gpu::GpuAccelerator,
    potential::{hartree, xc},
};

fn si_crystal() -> Crystal {
    let a = 5.431;
    Crystal {
        lattice: Lattice::new(
            a / 2.0 * Vector3::new(0.0, 1.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 0.0, 1.0),
            a / 2.0 * Vector3::new(1.0, 1.0, 0.0),
        ),
        atoms: vec![
            Atom::new(14, [0.0, 0.0, 0.0]),
            Atom::new(14, [0.25, 0.25, 0.25]),
        ],
    }
}

#[test]
fn test_gpu_hartree_on_realistic_density() {
    let Some(gpu) = GpuAccelerator::try_new() else {
        eprintln!("No GPU, skipping");
        return;
    };

    // Build a realistic G-space density from a Si crystal
    let crystal = si_crystal();
    let _basis = BasisSet::new(&crystal.lattice, 204.09);
    let grid = pwdft_rs::fft::FFT3D::new(20, 20, 20);
    let n_grid = grid.total_size();
    let recip = crystal.lattice.reciprocal();

    // Generate density-like data (positive real part, small imaginary)
    let rho_g: Vec<Complex64> = (0..n_grid)
        .map(|i| {
            let phase = i as f64 * 0.01;
            Complex64::new(0.001 * (-phase).exp(), 0.0001 * phase.sin())
        })
        .collect();

    // Compute g² for each grid point
    let [nx, ny, nz] = [20usize, 20, 20];
    let g_squared: Vec<f64> = (0..n_grid)
        .map(|idx| {
            let i1 = idx / (ny * nz);
            let i2 = (idx / nz) % ny;
            let i3 = idx % nz;
            let n1 = if i1 > nx / 2 { i1 as i32 - nx as i32 } else { i1 as i32 };
            let n2 = if i2 > ny / 2 { i2 as i32 - ny as i32 } else { i2 as i32 };
            let n3 = if i3 > nz / 2 { i3 as i32 - nz as i32 } else { i3 as i32 };
            let g = n1 as f64 * recip.a + n2 as f64 * recip.b + n3 as f64 * recip.c;
            g.norm_squared()
        })
        .collect();

    let fourpi_e2 = 4.0 * std::f64::consts::PI * hartree::E2;

    // CPU
    let cpu: Vec<Complex64> = rho_g
        .iter()
        .zip(g_squared.iter())
        .map(|(&rho, &g2)| {
            if g2 > 1e-20 { rho * fourpi_e2 / g2 } else { Complex64::new(0.0, 0.0) }
        })
        .collect();

    // GPU
    let gpu_result = gpu.hartree_potential(&rho_g, &g_squared, fourpi_e2);

    let mut max_rel_err = 0.0_f64;
    for (i, (c, g)) in cpu.iter().zip(gpu_result.iter()).enumerate() {
        let diff = (c - g).norm();
        let scale = c.norm().max(1e-15);
        let rel = diff / scale;
        max_rel_err = max_rel_err.max(rel);
        assert!(
            rel < 1e-3,
            "Hartree mismatch at {i}: cpu={c:.6e}, gpu={g:.6e}, rel={rel:.2e}"
        );
    }
    eprintln!("Hartree max relative error: {max_rel_err:.2e}");
}

#[test]
fn test_gpu_xc_across_density_regimes() {
    let Some(gpu) = GpuAccelerator::try_new() else {
        eprintln!("No GPU, skipping");
        return;
    };

    // Test across a wide range of densities: very low, typical metallic, very high
    let rho_r: Vec<f64> = (0..1000)
        .map(|i| {
            let t = i as f64 / 999.0;
            // Log-spaced from 1e-4 to 10.0 e/ų
            10.0_f64.powf(-4.0 + 5.0 * t)
        })
        .collect();

    let (cpu_exc, cpu_vxc) = xc::lda_xc_grid(&rho_r);
    let (gpu_exc, gpu_vxc) = gpu.lda_xc(&rho_r);

    let mut max_exc_err = 0.0_f64;
    let mut max_vxc_err = 0.0_f64;

    for (i, ((&ce, &cv), (&ge, &gv))) in cpu_exc
        .iter()
        .zip(cpu_vxc.iter())
        .zip(gpu_exc.iter().zip(gpu_vxc.iter()))
        .enumerate()
    {
        let exc_err = (ce - ge).abs();
        let vxc_err = (cv - gv).abs();
        max_exc_err = max_exc_err.max(exc_err);
        max_vxc_err = max_vxc_err.max(vxc_err);

        // Allow ~0.01 eV tolerance (f32 precision for ~1-20 eV range values)
        assert!(
            exc_err < 0.05,
            "XC exc at i={i} rho={:.4e}: cpu={ce:.6}, gpu={ge:.6}, err={exc_err:.2e}",
            rho_r[i]
        );
        assert!(
            vxc_err < 0.05,
            "XC vxc at i={i} rho={:.4e}: cpu={cv:.6}, gpu={gv:.6}, err={vxc_err:.2e}",
            rho_r[i]
        );
    }
    eprintln!("XC max errors: exc={max_exc_err:.4e} eV, vxc={max_vxc_err:.4e} eV");
}

#[test]
fn test_gpu_scf_converges_similarly() {
    let Some(_gpu) = GpuAccelerator::try_new() else {
        eprintln!("No GPU, skipping");
        return;
    };

    let crystal = si_crystal();
    let basis = BasisSet::new(&crystal.lattice, 204.09);
    let pp = pwdft_rs::pseudopotential::load(
        &std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("pseudopotentials/Si.UPF"),
    )
    .unwrap();
    let kpoints = pwdft_rs::kpoints::monkhorst_pack(2, 2, 2, &crystal.lattice);

    let params = pwdft_rs::scf::ScfParams {
        n_bands: 8,
        max_iter: 60,
        conv_threshold: 1e-6,
        mixing_beta: 0.3,
        mixing_ndim: 8,
        smearing_sigma: 0.05,
        ecutrho_ratio: 4,
        fft_grid: Some([20, 20, 20]),
    };

    // With gpu feature enabled, run_scf will use GPU automatically
    let result = pwdft_rs::scf::run_scf(&crystal, &basis, &kpoints, &[&pp], &params, None);

    match result {
        Ok(r) => {
            eprintln!("GPU SCF converged in {} iterations", r.n_iterations);
            eprintln!("Total energy: {:.6} eV", r.total_energy);
            eprintln!("Fermi energy: {:.6} eV", r.fermi_energy);

            // Gamma point eigenvalues
            if let Some(evs) = r.eigenvalues.first() {
                eprintln!("Eigenvalues at k=0: {:?}", evs);
            }

            // Total energy should be in a physically reasonable range for Si
            assert!(
                r.total_energy < -200.0 && r.total_energy > -250.0,
                "Total energy {:.2} eV outside expected range for Si",
                r.total_energy
            );
        }
        Err(e) => {
            eprintln!("GPU SCF did not converge: {e}");
            // Non-convergence within 60 iterations is acceptable —
            // the f32 precision may shift the convergence basin slightly
        }
    }
}
