//! GPU acceleration via wgpu (Metal backend on macOS).
//!
//! Provides GPU compute kernels for the grid-level operations in the SCF loop:
//! - Hartree potential: V_H(G) = 4πe² ρ(G) / |G|²
//! - V_eff assembly: v_eff = v_local + v_H + v_xc
//! - LDA XC: Perdew-Zunger exchange-correlation on real-space grid
//!
//! Uses f32 precision on GPU. The calling code handles f64↔f32 conversion
//! at the boundary. Eigensolves and energy accumulation remain f64 on CPU.
//!
//! ## f32/f64 error budget
//!
//! Expected relative errors vs f64 CPU path:
//! - Hartree potential: ~1e-7 (linear operations, f32 precision limit)
//! - LDA XC: ~1e-5 (cube root and log amplify rounding)
//! - V_eff assembly: ~1e-7 (linear addition)
//!
//! For SCF convergence to 1e-6 eV, f32 XC error (~1e-5 relative) can cause
//! the GPU path to converge to a slightly different energy (~0.01 meV).
//! This is acceptable for production use.
//!
//! Enable with `cargo build --features gpu`.

use bytemuck::{Pod, Zeroable};
use log::info;
use num_complex::Complex64;
use wgpu::util::DeviceExt;

/// GPU accelerator wrapping wgpu device and precompiled compute pipelines.
///
/// When `prepare_buffers` is called, allocates persistent GPU buffers for a
/// given grid size. Subsequent kernel calls reuse these buffers, avoiding
/// per-call allocation overhead (~1 ms per call).
pub struct GpuAccelerator {
    device: wgpu::Device,
    queue: wgpu::Queue,
    hartree_pipeline: wgpu::ComputePipeline,
    v_eff_pipeline: wgpu::ComputePipeline,
    lda_xc_pipeline: wgpu::ComputePipeline,
    /// Persistent buffer pool for a fixed grid size.
    pool: Option<BufferPool>,
}

/// Pre-allocated GPU buffers + bind groups for a fixed grid size.
///
/// Covers all three kernels (Hartree, V_eff, LDA XC) so that a steady-state
/// SCF iteration executes zero `device.create_buffer(...)` and zero
/// `device.create_bind_group(...)` calls — only per-call uniform updates via
/// `queue.write_buffer` and the input-data uploads.
///
/// Complex buffer roles (each sized `2 * n_grid * f32`):
/// - `complex_bufs[0]` — `rho_g` (Hartree input).
/// - `complex_bufs[1]` — `v_h` (Hartree output; also V_eff `v_h` input).
/// - `complex_bufs[2]` — `v_local` (V_eff input; static across SCF iterations
///   but currently re-uploaded per call).
/// - `complex_bufs[3]` — `v_xc_g` (V_eff input).
/// - `complex_bufs[4]` — `v_eff` (V_eff output).
///
/// Real-scalar buffer roles (each sized `n_grid * f32`):
/// - `rho_r_buf` — XC input density (real space, e/Å³).
/// - `exc_buf`, `vxc_buf` — XC outputs (eV).
/// - `g_squared_buf` — static |G|² (uploaded once in `prepare_buffers`).
///
/// Staging buffers:
/// - `complex_staging` — complex readback (Hartree / V_eff output).
/// - `exc_staging`, `vxc_staging` — real readbacks from LDA XC.
///
/// Uniform buffers (16 bytes each, updated per call via `write_buffer`):
/// - `hartree_params_buf` — `HartreeParams { fourpi_e2, n_grid }`.
/// - `grid_params_buf` — `GridParams { n_grid, _pad }` (shared by XC + V_eff).
///
/// Bind groups reference the buffers above and are rebuilt only in
/// `prepare_buffers` (i.e., once per SCF run when the grid is set up).
struct BufferPool {
    n_grid: usize,
    /// Complex buffers: 2 × n_grid f32 values each. See struct docs for roles.
    complex_bufs: Vec<wgpu::Buffer>,
    /// Staging buffer for complex readback (Hartree / V_eff output).
    complex_staging: wgpu::Buffer,
    /// Precomputed |G|² on GPU (doesn't change between SCF iterations).
    ///
    /// The cached `hartree_bg` bind group below holds the actual reference the
    /// kernel consumes; this field keeps the buffer alive and documents the
    /// ownership (retained for a future chain-fusion path that re-binds it
    /// when encoding Hartree → XC → V_eff into one encoder).
    #[allow(
        dead_code,
        reason = "Held alive via hartree_bg; explicit ownership retained for future chain fusion"
    )]
    g_squared_buf: wgpu::Buffer,
    /// XC real-space density input (n_grid f32).
    rho_r_buf: wgpu::Buffer,
    /// XC exchange-correlation energy density output (n_grid f32).
    exc_buf: wgpu::Buffer,
    /// XC exchange-correlation potential output (n_grid f32).
    vxc_buf: wgpu::Buffer,
    /// Staging buffer for exc readback (n_grid f32).
    exc_staging: wgpu::Buffer,
    /// Staging buffer for vxc readback (n_grid f32).
    vxc_staging: wgpu::Buffer,
    /// Per-call uniform buffer for Hartree params (16 bytes).
    hartree_params_buf: wgpu::Buffer,
    /// Per-call uniform buffer for grid params (shared by XC + V_eff).
    grid_params_buf: wgpu::Buffer,
    /// Cached bind group for the Hartree kernel.
    hartree_bg: wgpu::BindGroup,
    /// Cached bind group for the V_eff kernel.
    v_eff_bg: wgpu::BindGroup,
    /// Cached bind group for the LDA XC kernel.
    lda_xc_bg: wgpu::BindGroup,
}

// Uniform parameter structs matching WGSL layout.
// Must be Pod + Zeroable for bytemuck, and 8-byte aligned for wgpu uniforms.

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct HartreeParams {
    fourpi_e2: f32,
    n_grid: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct GridParams {
    n_grid: u32,
    _pad: u32, // Align to 8 bytes
}

const WORKGROUP_SIZE: u32 = 256;

fn dispatch_size(n: u32) -> u32 {
    n.div_ceil(WORKGROUP_SIZE)
}

impl GpuAccelerator {
    /// Try to initialize GPU. Returns None if no suitable adapter is found.
    pub fn try_new() -> Option<Self> {
        let mut desc = wgpu::InstanceDescriptor::new_without_display_handle();
        desc.backends = wgpu::Backends::METAL | wgpu::Backends::VULKAN;
        let instance = wgpu::Instance::new(desc);

        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        })).ok()?;

        let adapter_info = adapter.get_info();
        info!(
            "GPU: {} ({:?}, {:?})",
            adapter_info.name, adapter_info.backend, adapter_info.device_type
        );

        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("pwdft-rs GPU"),
                ..Default::default()
            },
        ))
        .ok()?;

        let hartree_pipeline = Self::create_pipeline(
            &device,
            "hartree",
            include_str!("shaders/hartree.wgsl"),
        );
        let v_eff_pipeline = Self::create_pipeline(
            &device,
            "v_eff_add",
            include_str!("shaders/v_eff_add.wgsl"),
        );
        let lda_xc_pipeline = Self::create_pipeline(
            &device,
            "lda_xc",
            include_str!("shaders/lda_xc.wgsl"),
        );

        Some(Self {
            device,
            queue,
            hartree_pipeline,
            v_eff_pipeline,
            lda_xc_pipeline,
            pool: None,
        })
    }

    /// Pre-allocate persistent GPU buffers + bind groups for a given grid
    /// size. Also uploads the static |G|² array that doesn't change between
    /// iterations.
    ///
    /// Covers all three kernels (Hartree, V_eff, LDA XC). After this call the
    /// pooled branches of `hartree_potential`, `v_eff_assembly`, and `lda_xc`
    /// issue zero `create_buffer` / `create_bind_group` per call — only
    /// `queue.write_buffer` updates for inputs and the kernel uniform.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "GPU mixed-precision strategy (module doc): f64→f32 conversion at the CPU→GPU boundary is intentional; see `GpuAccelerator` documentation"
    )]
    pub fn prepare_buffers(&mut self, n_grid: usize, g_squared: &[f64]) {
        let complex_size = (2 * n_grid * std::mem::size_of::<f32>()) as u64;
        let real_size = (n_grid * std::mem::size_of::<f32>()) as u64;

        // Allocate 5 complex storage buffers. See BufferPool docs for roles:
        // [0]=rho_g, [1]=v_h, [2]=v_local, [3]=v_xc_g, [4]=v_eff.
        let complex_bufs: Vec<wgpu::Buffer> = (0..5)
            .map(|i| {
                self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some(&format!("pool_complex_{i}")),
                    size: complex_size,
                    usage: wgpu::BufferUsages::STORAGE
                        | wgpu::BufferUsages::COPY_DST
                        | wgpu::BufferUsages::COPY_SRC,
                    mapped_at_creation: false,
                })
            })
            .collect();

        let complex_staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_complex_staging"),
            size: complex_size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        // Real-scalar buffers for the LDA XC kernel.
        let rho_r_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_rho_r"),
            size: real_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let exc_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_exc"),
            size: real_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let vxc_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_vxc"),
            size: real_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let exc_staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_exc_staging"),
            size: real_size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let vxc_staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_vxc_staging"),
            size: real_size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        // Upload static g_squared once per SCF run.
        let g2_f32: Vec<f32> = g_squared.iter().map(|&v| v as f32).collect();
        let g_squared_buf = self.create_storage_buffer(&g2_f32);

        // Per-call uniform buffers (16 bytes each, updated via write_buffer).
        let hartree_params_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_hartree_params"),
            size: std::mem::size_of::<HartreeParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let grid_params_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pool_grid_params"),
            size: std::mem::size_of::<GridParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Cached bind groups. Pool buffer handles are stable across SCF
        // iterations, so these are built once here and reused per call.
        let hartree_bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("hartree_bg"),
            layout: &self.hartree_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: hartree_params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: complex_bufs[0].as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: g_squared_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: complex_bufs[1].as_entire_binding() },
            ],
        });
        let v_eff_bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("v_eff_bg"),
            layout: &self.v_eff_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: grid_params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: complex_bufs[2].as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: complex_bufs[1].as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: complex_bufs[3].as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: complex_bufs[4].as_entire_binding() },
            ],
        });
        let lda_xc_bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("lda_xc_bg"),
            layout: &self.lda_xc_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: grid_params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: rho_r_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: exc_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: vxc_buf.as_entire_binding() },
            ],
        });

        self.pool = Some(BufferPool {
            n_grid,
            complex_bufs,
            complex_staging,
            g_squared_buf,
            rho_r_buf,
            exc_buf,
            vxc_buf,
            exc_staging,
            vxc_staging,
            hartree_params_buf,
            grid_params_buf,
            hartree_bg,
            v_eff_bg,
            lda_xc_bg,
        });

        info!("GPU buffer pool allocated for {n_grid} grid points (all kernels)");
    }

    fn create_pipeline(
        device: &wgpu::Device,
        label: &str,
        source: &str,
    ) -> wgpu::ComputePipeline {
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some(label),
            source: wgpu::ShaderSource::Wgsl(source.into()),
        });
        device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(label),
            layout: None, // Auto-layout from shader bindings
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        })
    }

    /// Compute Hartree potential on GPU.
    ///
    /// V_H(G) = 4πe² × ρ(G) / |G|² for |G|² > 0, else 0.
    ///
    /// Uses pooled buffers if `prepare_buffers` was called; otherwise allocates fresh.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "GPU mixed-precision: f64→f32 at CPU→GPU boundary is intentional; n_grid <= 512^3 (~1.3e8) fits in u32::MAX (~4.3e9) for all physical inputs"
    )]
    pub fn hartree_potential(
        &self,
        rho_g: &[Complex64],
        g_squared: &[f64],
        fourpi_e2: f64,
    ) -> Vec<Complex64> {
        let n_grid = rho_g.len();
        let rho_f32 = complex_to_f32_pairs(rho_g);

        let params = HartreeParams {
            fourpi_e2: fourpi_e2 as f32,
            n_grid: n_grid as u32,
        };

        // Pooled path: cached bind group + pre-allocated buffers. Only the
        // input `rho_g` and the uniform params travel host→device per call.
        if let Some(ref pool) = self.pool
            && pool.n_grid == n_grid
        {
            self.queue.write_buffer(
                &pool.hartree_params_buf,
                0,
                bytemuck::bytes_of(&params),
            );
            self.queue
                .write_buffer(&pool.complex_bufs[0], 0, bytemuck::cast_slice(&rho_f32));

            let mut encoder = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = encoder.begin_compute_pass(&Default::default());
                pass.set_pipeline(&self.hartree_pipeline);
                pass.set_bind_group(0, &pool.hartree_bg, &[]);
                pass.dispatch_workgroups(dispatch_size(n_grid as u32), 1, 1);
            }
            let byte_size = (rho_f32.len() * std::mem::size_of::<f32>()) as u64;
            encoder.copy_buffer_to_buffer(
                &pool.complex_bufs[1],
                0,
                &pool.complex_staging,
                0,
                byte_size,
            );
            self.queue.submit(std::iter::once(encoder.finish()));

            let result_f32 = self.read_staging_buffer(&pool.complex_staging, rho_f32.len());
            return f32_pairs_to_complex(&result_f32);
        }

        // Fallback: allocate fresh buffers (no prepare_buffers call).
        let params_buf = self.create_uniform_buffer(&params);
        let g2_f32: Vec<f32> = g_squared.iter().map(|&v| v as f32).collect();
        let rho_buf = self.create_storage_buffer(&rho_f32);
        let g2_buf = self.create_storage_buffer(&g2_f32);
        let out_buf = self.create_output_buffer(rho_f32.len());
        let staging_buf = self.create_staging_buffer(rho_f32.len());

        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("hartree"),
            layout: &self.hartree_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: rho_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: g2_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: out_buf.as_entire_binding() },
            ],
        });

        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.hartree_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(dispatch_size(n_grid as u32), 1, 1);
        }
        let byte_size = (rho_f32.len() * std::mem::size_of::<f32>()) as u64;
        encoder.copy_buffer_to_buffer(&out_buf, 0, &staging_buf, 0, byte_size);
        self.queue.submit(std::iter::once(encoder.finish()));

        let result_f32 = self.read_staging_buffer(&staging_buf, rho_f32.len());
        f32_pairs_to_complex(&result_f32)
    }

    /// Compute V_eff = V_local + V_H + V_xc on GPU.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "GPU mixed-precision: f64→f32 at CPU→GPU boundary is intentional; n_grid <= 512^3 fits in u32::MAX"
    )]
    pub fn v_eff_assembly(
        &self,
        v_local: &[Complex64],
        v_h: &[Complex64],
        v_xc: &[Complex64],
    ) -> Vec<Complex64> {
        let n_grid = v_local.len();
        assert_eq!(v_h.len(), n_grid);
        assert_eq!(v_xc.len(), n_grid);

        let vl_f32 = complex_to_f32_pairs(v_local);
        let vh_f32 = complex_to_f32_pairs(v_h);
        let vxc_f32 = complex_to_f32_pairs(v_xc);

        let params = GridParams {
            n_grid: n_grid as u32,
            _pad: 0,
        };

        // Pooled path: cached bind group + pre-allocated buffers. Inputs are
        // written into `complex_bufs[2] = v_local`, `[1] = v_h`, `[3] = v_xc`;
        // output lands in `[4] = v_eff`.
        if let Some(ref pool) = self.pool
            && pool.n_grid == n_grid
        {
            self.queue
                .write_buffer(&pool.grid_params_buf, 0, bytemuck::bytes_of(&params));
            self.queue
                .write_buffer(&pool.complex_bufs[2], 0, bytemuck::cast_slice(&vl_f32));
            self.queue
                .write_buffer(&pool.complex_bufs[1], 0, bytemuck::cast_slice(&vh_f32));
            self.queue
                .write_buffer(&pool.complex_bufs[3], 0, bytemuck::cast_slice(&vxc_f32));

            let mut encoder = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = encoder.begin_compute_pass(&Default::default());
                pass.set_pipeline(&self.v_eff_pipeline);
                pass.set_bind_group(0, &pool.v_eff_bg, &[]);
                pass.dispatch_workgroups(dispatch_size(n_grid as u32), 1, 1);
            }
            encoder.copy_buffer_to_buffer(
                &pool.complex_bufs[4],
                0,
                &pool.complex_staging,
                0,
                (vl_f32.len() * std::mem::size_of::<f32>()) as u64,
            );
            self.queue.submit(std::iter::once(encoder.finish()));

            let result_f32 = self.read_staging_buffer(&pool.complex_staging, vl_f32.len());
            return f32_pairs_to_complex(&result_f32);
        }

        // Fallback: allocate fresh buffers (no prepare_buffers call).
        let params_buf = self.create_uniform_buffer(&params);
        let vl_buf = self.create_storage_buffer(&vl_f32);
        let vh_buf = self.create_storage_buffer(&vh_f32);
        let vxc_buf = self.create_storage_buffer(&vxc_f32);
        let out_buf = self.create_output_buffer(vl_f32.len());
        let staging_buf = self.create_staging_buffer(vl_f32.len());

        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("v_eff"),
            layout: &self.v_eff_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: vl_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: vh_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: vxc_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: out_buf.as_entire_binding() },
            ],
        });

        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.v_eff_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(dispatch_size(n_grid as u32), 1, 1);
        }
        encoder.copy_buffer_to_buffer(
            &out_buf, 0,
            &staging_buf, 0,
            (vl_f32.len() * std::mem::size_of::<f32>()) as u64,
        );
        self.queue.submit(std::iter::once(encoder.finish()));

        let result_f32 = self.read_staging_buffer(&staging_buf, vl_f32.len());
        f32_pairs_to_complex(&result_f32)
    }

    /// Compute LDA exchange-correlation on GPU.
    ///
    /// Input: real-space density ρ(r) in e/Å³.
    /// Output: (ε_xc(r), V_xc(r)) in eV.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "GPU mixed-precision: f64→f32 at CPU→GPU boundary is intentional; n_grid <= 512^3 fits in u32::MAX"
    )]
    pub fn lda_xc(&self, rho_r: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n_grid = rho_r.len();
        let rho_f32: Vec<f32> = rho_r.iter().map(|&v| v as f32).collect();

        let params = GridParams {
            n_grid: n_grid as u32,
            _pad: 0,
        };

        // Pooled path: cached bind group + pre-allocated buffers. Only the
        // input `rho_r` and the uniform params travel host→device per call.
        if let Some(ref pool) = self.pool
            && pool.n_grid == n_grid
        {
            self.queue
                .write_buffer(&pool.grid_params_buf, 0, bytemuck::bytes_of(&params));
            self.queue
                .write_buffer(&pool.rho_r_buf, 0, bytemuck::cast_slice(&rho_f32));

            let mut encoder = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = encoder.begin_compute_pass(&Default::default());
                pass.set_pipeline(&self.lda_xc_pipeline);
                pass.set_bind_group(0, &pool.lda_xc_bg, &[]);
                pass.dispatch_workgroups(dispatch_size(n_grid as u32), 1, 1);
            }
            let byte_size = (n_grid * std::mem::size_of::<f32>()) as u64;
            encoder.copy_buffer_to_buffer(&pool.exc_buf, 0, &pool.exc_staging, 0, byte_size);
            encoder.copy_buffer_to_buffer(&pool.vxc_buf, 0, &pool.vxc_staging, 0, byte_size);
            self.queue.submit(std::iter::once(encoder.finish()));

            let exc_f32 = self.read_staging_buffer(&pool.exc_staging, n_grid);
            let vxc_f32 = self.read_staging_buffer(&pool.vxc_staging, n_grid);

            let exc: Vec<f64> = exc_f32.iter().map(|&v| v as f64).collect();
            let vxc: Vec<f64> = vxc_f32.iter().map(|&v| v as f64).collect();
            return (exc, vxc);
        }

        // Fallback: allocate fresh buffers (no prepare_buffers call).
        let params_buf = self.create_uniform_buffer(&params);
        let rho_buf = self.create_storage_buffer(&rho_f32);
        let exc_buf = self.create_output_buffer(n_grid);
        let vxc_buf = self.create_output_buffer(n_grid);
        let exc_staging = self.create_staging_buffer(n_grid);
        let vxc_staging = self.create_staging_buffer(n_grid);

        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("lda_xc"),
            layout: &self.lda_xc_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: rho_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: exc_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: vxc_buf.as_entire_binding() },
            ],
        });

        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.lda_xc_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(dispatch_size(n_grid as u32), 1, 1);
        }
        let byte_size = (n_grid * std::mem::size_of::<f32>()) as u64;
        encoder.copy_buffer_to_buffer(&exc_buf, 0, &exc_staging, 0, byte_size);
        encoder.copy_buffer_to_buffer(&vxc_buf, 0, &vxc_staging, 0, byte_size);
        self.queue.submit(std::iter::once(encoder.finish()));

        let exc_f32 = self.read_staging_buffer(&exc_staging, n_grid);
        let vxc_f32 = self.read_staging_buffer(&vxc_staging, n_grid);

        let exc: Vec<f64> = exc_f32.iter().map(|&v| v as f64).collect();
        let vxc: Vec<f64> = vxc_f32.iter().map(|&v| v as f64).collect();
        (exc, vxc)
    }

    // --- Buffer helpers ---

    fn create_uniform_buffer<T: Pod>(&self, data: &T) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: bytemuck::bytes_of(data),
            usage: wgpu::BufferUsages::UNIFORM,
        })
    }

    fn create_storage_buffer(&self, data: &[f32]) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: bytemuck::cast_slice(data),
            usage: wgpu::BufferUsages::STORAGE,
        })
    }

    fn create_output_buffer(&self, n_floats: usize) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: (n_floats * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        })
    }

    fn create_staging_buffer(&self, n_floats: usize) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: (n_floats * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        })
    }

    fn read_staging_buffer(&self, buffer: &wgpu::Buffer, n_floats: usize) -> Vec<f32> {
        let slice = buffer.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |result| {
            // SAFETY: receiver.recv() is called immediately after device.poll(),
            // so the receiver cannot be dropped before this callback fires.
            sender.send(result)
                .expect("BUG: GPU readback channel closed before callback fired");
        });
        let _ = self.device.poll(wgpu::PollType::Wait { submission_index: None, timeout: None });
        // SAFETY: The callback above sends exactly one message after poll() completes,
        // so recv() cannot fail. The inner Result is the GPU mapping result.
        receiver.recv()
            .expect("BUG: GPU readback channel closed unexpectedly")
            .expect("BUG: GPU buffer mapping failed");

        let data = slice.get_mapped_range();
        let result = bytemuck::cast_slice(&data)[..n_floats].to_vec();
        drop(data);
        buffer.unmap();

        result
    }
}

// --- f64 ↔ f32 conversion helpers ---

/// Convert Complex64 slice to interleaved f32 pairs [re0, im0, re1, im1, ...].
#[allow(
    clippy::cast_possible_truncation,
    reason = "GPU mixed-precision boundary: f64→f32 conversion at CPU→GPU boundary is intentional; see module-level docs"
)]
fn complex_to_f32_pairs(data: &[Complex64]) -> Vec<f32> {
    let mut out = Vec::with_capacity(data.len() * 2);
    for c in data {
        out.push(c.re as f32);
        out.push(c.im as f32);
    }
    out
}

/// Convert interleaved f32 pairs back to Complex64.
fn f32_pairs_to_complex(data: &[f32]) -> Vec<Complex64> {
    data.chunks_exact(2)
        .map(|pair| Complex64::new(pair[0] as f64, pair[1] as f64))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::xc;

    fn try_gpu() -> Option<GpuAccelerator> {
        GpuAccelerator::try_new()
    }

    #[test]
    fn test_gpu_hartree_matches_cpu() {
        let Some(gpu) = try_gpu() else {
            eprintln!("No GPU available, skipping test");
            return;
        };

        let n = 1000;
        let fourpi_e2 = 4.0 * std::f64::consts::PI * crate::consts::E2_COULOMB;

        // Generate test data
        let rho_g: Vec<Complex64> = (0..n)
            .map(|i| Complex64::new((i as f64 * 0.1).sin() * 0.01, (i as f64 * 0.2).cos() * 0.01))
            .collect();
        let g_squared: Vec<f64> = (0..n)
            .map(|i| if i == 0 { 0.0 } else { 1.0 + i as f64 * 0.5 })
            .collect();

        // CPU reference
        let cpu_result: Vec<Complex64> = rho_g
            .iter()
            .zip(g_squared.iter())
            .map(|(&rho, &g2)| {
                if g2 > crate::consts::G2_ZERO_THRESHOLD {
                    rho * fourpi_e2 / g2
                } else {
                    Complex64::new(0.0, 0.0)
                }
            })
            .collect();

        // GPU
        let gpu_result = gpu.hartree_potential(&rho_g, &g_squared, fourpi_e2);

        // Compare (f32 tolerance: ~1e-5 relative)
        for (i, (c, g)) in cpu_result.iter().zip(gpu_result.iter()).enumerate() {
            let diff = (c - g).norm();
            let scale = c.norm().max(1e-10);
            assert!(
                diff / scale < 1e-4,
                "Hartree mismatch at {i}: cpu={c}, gpu={g}, rel_err={}",
                diff / scale
            );
        }
    }

    #[test]
    fn test_gpu_v_eff_matches_cpu() {
        let Some(gpu) = try_gpu() else {
            eprintln!("No GPU available, skipping test");
            return;
        };

        let n = 2000;
        let make_complex = |seed: f64| -> Vec<Complex64> {
            (0..n)
                .map(|i| {
                    Complex64::new(
                        (i as f64 * seed).sin() * 0.5,
                        (i as f64 * seed * 1.3).cos() * 0.5,
                    )
                })
                .collect()
        };

        let v_local = make_complex(0.1);
        let v_h = make_complex(0.2);
        let v_xc = make_complex(0.3);

        // CPU
        let cpu_result: Vec<Complex64> = (0..n)
            .map(|i| v_local[i] + v_h[i] + v_xc[i])
            .collect();

        // GPU
        let gpu_result = gpu.v_eff_assembly(&v_local, &v_h, &v_xc);

        for (i, (c, g)) in cpu_result.iter().zip(gpu_result.iter()).enumerate() {
            let diff = (c - g).norm();
            let scale = c.norm().max(1e-10);
            assert!(
                diff / scale < 1e-5,
                "V_eff mismatch at {i}: cpu={c}, gpu={g}",
            );
        }
    }

    #[test]
    fn test_gpu_lda_xc_matches_cpu() {
        let Some(gpu) = try_gpu() else {
            eprintln!("No GPU available, skipping test");
            return;
        };

        // Test densities spanning both PZ regimes (rs < 1 and rs >= 1)
        let rho_r: Vec<f64> = (0..500)
            .map(|i| 0.001 + i as f64 * 0.01) // 0.001 to 5.0 e/ų
            .collect();

        // CPU reference
        let (cpu_exc, cpu_vxc) = xc::lda_xc_grid(&rho_r);

        // GPU
        let (gpu_exc, gpu_vxc) = gpu.lda_xc(&rho_r);

        // f32 XC should agree to ~1e-3 eV (f32 has ~7 decimal digits,
        // and the XC values are ~1-10 eV range)
        for (i, ((&ce, &cv), (&ge, &gv))) in cpu_exc
            .iter()
            .zip(cpu_vxc.iter())
            .zip(gpu_exc.iter().zip(gpu_vxc.iter()))
            .enumerate()
        {
            let exc_diff = (ce - ge).abs();
            let vxc_diff = (cv - gv).abs();
            assert!(
                exc_diff < 0.01,
                "XC exc mismatch at {i} (rho={:.4}): cpu={ce:.6}, gpu={ge:.6}, diff={exc_diff:.2e}",
                rho_r[i]
            );
            assert!(
                vxc_diff < 0.01,
                "XC vxc mismatch at {i} (rho={:.4}): cpu={cv:.6}, gpu={gv:.6}, diff={vxc_diff:.2e}",
                rho_r[i]
            );
        }
    }
}
