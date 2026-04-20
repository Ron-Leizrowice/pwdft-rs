// LDA exchange-correlation: Perdew-Zunger parametrization.
// Input:  rho_r (real-space density, f32, in e/Å³)
// Output: exc_r (energy density, eV), vxc_r (potential, eV)
//
// References:
//   Exchange: Slater, Phys. Rev. 81, 385 (1951)
//   Correlation: Perdew & Zunger, Phys. Rev. B 23, 5048 (1981), Table I
//
// Runs in f32 for GPU throughput. Expected relative error vs f64 CPU: ~1e-5
// (dominated by cube root and log rounding in f32).

struct Params {
    n_grid: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> rho_r: array<f32>;
@group(0) @binding(2) var<storage, read_write> exc_r: array<f32>;
@group(0) @binding(3) var<storage, read_write> vxc_r: array<f32>;

const PI: f32 = 3.14159265358979323846;
const HA_TO_EV: f32 = 27.211386;
const BOHR3: f32 = 0.1481847; // 0.529177^3 (7 significant digits, f32 limit)

// Perdew-Zunger correlation parameters (unpolarized, rs >= 1)
const PZ_GAMMA: f32 = -0.1423;
const PZ_BETA1: f32 = 1.0529;
const PZ_BETA2: f32 = 0.3334;

// Perdew-Zunger correlation parameters (unpolarized, rs < 1)
const PZ_A: f32 = 0.0311;
const PZ_B: f32 = -0.048;
const PZ_C: f32 = 0.0020;
const PZ_D: f32 = -0.0116;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if (idx >= params.n_grid) {
        return;
    }

    let rho = rho_r[idx];
    // Must match RHO_FLOOR in src/consts.rs
    if (rho < 1e-20) {
        exc_r[idx] = 0.0;
        vxc_r[idx] = 0.0;
        return;
    }

    // Convert to Bohr units
    let rho_bohr = rho * BOHR3;

    // Wigner-Seitz radius
    let rs = pow(3.0 / (4.0 * PI * rho_bohr), 1.0 / 3.0);

    // Slater exchange: ε_x = -(3/4)(3ρ/π)^{1/3} in Hartree
    let cbrt_arg = pow(3.0 * rho_bohr / PI, 1.0 / 3.0);
    let ex_ha = -0.75 * cbrt_arg;
    let vx_ha = (4.0 / 3.0) * ex_ha;

    // Perdew-Zunger correlation
    var ec_ha: f32;
    var vc_ha: f32;

    if (rs >= 1.0) {
        let sqrt_rs = sqrt(rs);
        let denom = 1.0 + PZ_BETA1 * sqrt_rs + PZ_BETA2 * rs;
        ec_ha = PZ_GAMMA / denom;
        let d_ec = -PZ_GAMMA * (PZ_BETA1 / (2.0 * sqrt_rs) + PZ_BETA2) / (denom * denom);
        vc_ha = ec_ha - rs / 3.0 * d_ec;
    } else {
        let ln_rs = log(rs);
        ec_ha = PZ_A * ln_rs + PZ_B + PZ_C * rs * ln_rs + PZ_D * rs;
        let d_ec = PZ_A / rs + PZ_C * (ln_rs + 1.0) + PZ_D;
        vc_ha = ec_ha - rs / 3.0 * d_ec;
    }

    exc_r[idx] = (ex_ha + ec_ha) * HA_TO_EV;
    vxc_r[idx] = (vx_ha + vc_ha) * HA_TO_EV;
}
