// Hartree potential: V_H(G) = 4πe² × ρ(G) / |G|²
// Complex numbers stored as consecutive f32 pairs: [re0, im0, re1, im1, ...]

struct Params {
    fourpi_e2: f32,
    n_grid: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> rho_g: array<f32>;
@group(0) @binding(2) var<storage, read> g_squared: array<f32>;
@group(0) @binding(3) var<storage, read_write> v_h: array<f32>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if (idx >= params.n_grid) {
        return;
    }

    let g2 = g_squared[idx];
    let i_re = 2u * idx;
    let i_im = i_re + 1u;

    if (g2 > 1e-20) {
        let factor = params.fourpi_e2 / g2;
        v_h[i_re] = rho_g[i_re] * factor;
        v_h[i_im] = rho_g[i_im] * factor;
    } else {
        v_h[i_re] = 0.0;
        v_h[i_im] = 0.0;
    }
}
