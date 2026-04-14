// V_eff assembly: v_eff(G) = v_local(G) + v_H(G) + v_xc(G)
// All arrays are complex (f32 pairs).

struct Params {
    n_grid: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> v_local: array<f32>;
@group(0) @binding(2) var<storage, read> v_h: array<f32>;
@group(0) @binding(3) var<storage, read> v_xc: array<f32>;
@group(0) @binding(4) var<storage, read_write> v_eff: array<f32>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if (idx >= params.n_grid) {
        return;
    }

    let i_re = 2u * idx;
    let i_im = i_re + 1u;

    v_eff[i_re] = v_local[i_re] + v_h[i_re] + v_xc[i_re];
    v_eff[i_im] = v_local[i_im] + v_h[i_im] + v_xc[i_im];
}
