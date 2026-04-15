# Pseudopotential approaches for a new plane-wave DFT code

**PAW remains the accuracy gold standard, but starting with ONCV norm-conserving pseudopotentials is the strategically optimal path for a new Rust-based plane-wave code.** Modern ONCV libraries (PseudoDojo, SG15, SPMS) now achieve accuracy within ~0.3–1.0 meV/atom of all-electron results while requiring dramatically simpler implementation than PAW — no augmentation charges, no generalized eigenvalue problem, no one-center corrections. A phased approach — ONCV first via UPF2 parsing, PAW second via PAW-XML — maximizes early capability while building toward VASP-competitive accuracy. For actinide systems like UO₂, the 2024 PseudoDojo actinide extension by Tantardini et al. provides the first comprehensive open ONCV library covering Z = 87–118, making norm-conserving calculations on uranium feasible without proprietary VASP potentials.

---

## The three formalisms and why PAW dominates — but shouldn't come first

Three pseudopotential formalisms compete in modern plane-wave DFT. **Norm-conserving pseudopotentials (NCPP)** replace core electrons with a smooth effective potential while preserving the norm of pseudo-wavefunctions beyond a cutoff radius. The original Hamann-Schlüter-Chiang formulation (Phys. Rev. Lett. 43, 1494, 1979) was refined by Troullier and Martins (Phys. Rev. B 43, 1993, 1991) for smoothness, and then fundamentally advanced by Hamann's Optimized Norm-Conserving Vanderbilt (ONCV) method (Phys. Rev. B 88, 085117, 2013), which uses two projectors per angular momentum channel to achieve softness competitive with ultrasoft potentials while retaining the simplicity of norm conservation. ONCV potentials typically require **60–100 Ry** wavefunction cutoffs — roughly 2× higher than PAW/USPP, but manageable on modern hardware and dramatically lower than classical NCPPs.

**Ultrasoft pseudopotentials (USPP)**, introduced by Vanderbilt (Phys. Rev. B 41, 7892, 1990), relax the norm-conservation constraint to achieve much softer potentials (**30–50 Ry** cutoffs), at the cost of introducing augmentation charges and a generalized eigenvalue problem. The GBRV library (Garrity et al., Comput. Mater. Sci. 81, 446, 2014) demonstrated RMS lattice constant errors below 0.15% versus all-electron WIEN2k. However, USPP development has largely stalled — the GBRV library (v1.5) is the last major release, and USPPs cannot reconstruct true all-electron wavefunctions near nuclei.

**The Projector Augmented Wave (PAW) method** by Blöchl (Phys. Rev. B 50, 17953, 1994) defines a linear transformation from smooth pseudo-wavefunctions to true all-electron wavefunctions via partial waves and projectors within augmentation spheres. Kresse and Joubert (Phys. Rev. B 59, 1758, 1999) showed that USPP is formally a linearized approximation to PAW, and their VASP implementation made PAW the de facto standard. PAW achieves **all-electron accuracy within the frozen-core approximation** at USPP-like cutoffs, provides access to properties requiring full wavefunctions (NMR, EPR, core-level spectroscopy), and is systematically improvable via additional projectors.

| Feature | NCPP (ONCV) | USPP | PAW |
|---|---|---|---|
| Typical wavefunction cutoff | 60–100 Ry | 30–50 Ry | 30–50 Ry |
| Eigenvalue problem | Standard | Generalized | Generalized |
| Augmentation charges | None | Required | Required |
| All-electron reconstruction | No | No | Yes |
| Implementation complexity | Lowest | Medium | Highest |
| Active development | Strong | Limited | Strong |

**For a new code, the recommendation is unambiguous: implement ONCV first.** The implementation requires only local potential application, Kleinman-Bylander nonlocal projectors, and the D_ij coefficient matrix — no overlap operator, no augmentation machinery. PAW adds the overlap operator S, augmentation charges Q_ij, compensation charges with shape functions, the three-term energy decomposition (Ẽ + E¹ − Ẽ¹), and modified eigensolvers for the generalized eigenvalue problem. This roughly triples implementation effort. Modern ONCV potentials are accurate enough for most applications, and NCPP is *required* for several advanced methods including GW, Wannier functions, TDDFT, and certain DFPT implementations.

---

## Libraries and benchmarks: the open-source ecosystem has matured

The pseudopotential library landscape has consolidated around a few high-quality, well-benchmarked collections. The primary verification metric is the **Δ-gauge** (Lejaeghere et al., Science 351, aad3000, 2016), which measures the RMS difference between equation-of-state curves from pseudopotential and all-electron calculations across 71 elemental crystals, expressed in meV/atom. Modern pseudopotential codes achieve **Δ values of 0.3–1.0 meV/atom**, comparable to differences between all-electron codes themselves.

**PseudoDojo** (van Setten et al., Comput. Phys. Commun. 226, 39, 2018) is the most comprehensive open library, providing ONCV pseudopotentials for **85 elements** in two stringency levels — "standard" (softer, lower cutoffs) and "stringent" (harder, higher accuracy) — across PBE, PBEsol, and LDA functionals. Both scalar-relativistic and fully-relativistic versions are available. Formats include psp8, UPF2, and PSML 1.1. The library underwent seven batteries of automated tests via ABINIT, and the 2024 actinide extension (Tantardini et al., Comput. Phys. Commun. 295, 108986) added 34 pseudopotentials for elements Z = 87–118, generated with ONCVPSP v4.0.1 and validated against all-electron ZORA calculations. PseudoDojo is the **recommended primary target** for a new code.

**SSSP (Standard Solid-State Pseudopotentials)** from Materials Cloud (Prandini et al., npj Comput. Mater. 4, 72, 2018) takes a meta-library approach, selecting the best pseudopotential for each element from PseudoDojo, GBRV, PSlibrary, and SG15 based on Δ-factor testing, phonon convergence, band structure convergence, and cohesive energy convergence. Two protocols exist: **SSSP Efficiency** (optimized for throughput) and **SSSP Precision** (highest accuracy, scoring best among open-source libraries on the Δ-factor test). SSSP mixes NC, USPP, and PAW potentials — a code supporting only NC would access a subset.

**SG15** (Schlipf and Gygi, Comput. Phys. Commun. 196, 36, 2015) provides ONCV potentials for elements up to Z = 83 (excluding lanthanides), automatically optimized via a Nelder-Mead algorithm balancing accuracy and softness. Lattice constant deviations of ~0.1% versus FLAPW are typical at ~60 Ry cutoff. **GBRV** (v1.5) remains the standard ultrasoft library but excludes f-block elements. **JTH PAW tables** (Jollet, Torrent, Holzwarth, Comput. Phys. Commun. 185, 1246, 2014) provide PAW-XML datasets for ABINIT/GPAW, currently at v2.0 with coverage of 71+ elements including rare earths.

The **SPMS library** (Shojaei et al., Comput. Phys. Commun. 284, 108594, 2023) deserves special attention: it uses multi-objective evolutionary optimization to generate ONCV potentials that are **~35% softer** than PseudoDojo at comparable accuracy (average cutoff 18.7 Ha versus 29.1 Ha), yielding 2–5× speedups. Available at github.com/SPARC-X/SPMS-psps for 69 elements.

The landmark **Bosoni et al.** study (Nat. Rev. Phys. 6, 45, 2024) expanded verification to Z = 1–96 with 960 equations of state across 10 prototypical cubic compounds per element, establishing the definitive benchmark framework using AiiDA reproducible workflows.

---

## Relativistic treatments for heavy elements

Relativistic effects are incorporated at three levels. **Scalar-relativistic** treatment, standard in all modern pseudopotential generators, uses the Koelling-Harmon approach (J. Phys. C 10, 3107, 1977) to include mass-velocity and Darwin terms while omitting spin-orbit coupling. This captures orbital contraction/expansion and energy level shifts — effects that can change bond lengths by several percent for heavy elements. All VASP PAW potentials, PseudoDojo pseudopotentials, and JTH datasets are generated scalar-relativistically by default.

**Fully-relativistic pseudopotentials** solve the four-component Dirac equation for the atom, yielding j-dependent projectors split into j = l+1/2 and j = l−1/2 channels (Theurich and Hill, Phys. Rev. B 64, 073106, 2001). In plane-wave calculations, these are applied via a two-component spinor formalism — the small Dirac component has negligible density in the valence region. SOC enters as ΔW_l^SO(r) = W_{l,j=l+1/2}(r) − W_{l,j=l−1/2}(r). PseudoDojo provides fully-relativistic versions for all 85+ elements; SG15 includes FR potentials generated by Scherpelz. For topological insulators, heavy-element semiconductors (Pb halide perovskites: SOC reduces band gaps by ~1 eV), and all actinide/lanthanide systems, **FR pseudopotentials are essential**.

Alternative relativistic approximations in pseudopotential generation include **ZORA** (Zeroth-Order Regular Approximation), which avoids Pauli-expansion divergences and is used in CASTEP and AMS/BAND, and the **exact Dirac equation** solution used by ONCVPSP and most modern generators. For elements with Z > 50, the choice between these approaches primarily affects SOC splittings and band gaps.

**Semicore states** — inner-shell electrons with significant spatial overlap with valence orbitals — must be included explicitly for transition metals (3s3p for 3d metals), lanthanides (5s5p, sometimes 4f), and actinides (6s6p always, sometimes 5d). VASP labels these as _pv (p-semicore), _sv (s+p semicore), with different ENMAX cutoffs. Including semicore states typically increases cutoffs by 50–200 eV but is mandatory for DFT+U calculations, high-pressure studies, and any system where core-valence hybridization occurs. For uranium, the standard 14-electron configuration (6s²6p⁶5f³6d¹7s²) is universally adopted.

---

## Actinide pseudopotentials and the UO₂ challenge

Uranium dioxide is a prototypical strongly-correlated insulator where standard DFT fails catastrophically — LDA/GGA incorrectly predicts it as metallic. The U⁴⁺ ion has a 5f² configuration with strongly localized electrons, making it a Mott-Hubbard insulator with an experimental band gap of ~2 eV. Three beyond-DFT approaches address this, all critically dependent on pseudopotential quality.

**DFT+U** is the most widely used approach. The Hubbard correction is applied to 5f orbital projectors defined by the pseudopotential, making results directly sensitive to PP choice. Published UO₂ studies use U values of **3.5–4.5 eV** for U-5f, but the optimal value depends on the pseudopotential: Allen and Watson (J. Phys. Chem. C 126, 12247, 2022) recommend U = 4.0 eV with PBESol PAW in VASP; Dorado et al. (Phys. Rev. B 79, 235125, 2009) used U = 4.50 eV, J = 0.54 eV; a 2025 QE study found U = 2.45 eV via self-consistent linear response gives the correct lattice parameter (5.47 Å). **Occupation matrix control (OMC)** is essential — Dorado et al. demonstrated that without it, DFT+U converges to metastable states spanning 3.45 eV in energy and ranging from metallic to 2.8 eV band gap across 21 possible f-orbital configurations. Alternative approaches include U-ramping (Krack, CP2K) and controlled symmetry reduction.

**Hybrid functionals** (HSE06 with 25% short-range HF exchange, ω = 0.11 Bohr⁻¹) give UO₂ band gaps closer to experiment without empirical U parameters. VASP explicitly warns against using soft (_s) potentials for hybrid calculations. Computational cost scales as N_k² × N_bands² × N_G log(N_G), making these calculations expensive but increasingly routine.

**DFT+DMFT** requires localized projector functions for the correlated 5f subspace. PAW naturally provides atomic-like projectors suitable for this purpose. ABINIT has a built-in DFT+DMFT implementation using projected Wannier functions (Amadon, J. Phys.: Condens. Matter 20, 235210, 2008), while VASP interfaces with TRIQS/DFTTools via LOCPROJ/PLOVASP. The PAW formalism is strongly preferred for DMFT workflows.

For pseudopotential availability, the 2024 PseudoDojo actinide extension by Tantardini et al. is the most significant recent development. It provides **ONCV pseudopotentials for all actinides (Ac–Lr) plus super-heavy elements (Rf–Og)**, generated with ONCVPSP v4.0.1 in scalar-relativistic and fully-relativistic versions across PBE, PBEsol, and LDA. These were validated against all-electron ZORA calculations using Δ-gauge and Δ₁-gauge descriptors. The paper explicitly notes that "the chemical bonding of actinides cannot be described if 5f-electrons are frozen in the core." VASP PAW potentials cover Ac–Pu but are proprietary. The JTH v2.0 PAW table extends through Pu. CP2K offers GTH potentials for the full actinide series in medium-core and large-core variants (Cantu et al., ~2021).

---

## File formats: UPF2 first, PSP8 second, PAW-XML for the PAW phase

Five pseudopotential file formats matter. **UPF (Unified Pseudopotential Format)** is the highest-priority target. Version 2.0.1 is nominally XML but uses non-standard conventions (enumerated tags like `<PP_BETA.1>`, Fortran-style booleans) that break standard XML parsers. It stores local potential, Kleinman-Bylander projectors with D_ij matrix, augmentation charges (for USPP/PAW), pseudo-wavefunctions, atomic charge density, and optional sections for PAW (all-electron partial waves, augmentation multipoles), spin-orbit data, and GIPAW reconstruction. Units are atomic Rydberg. The format supports NC, USPP, and PAW in a single specification. PseudoDojo, SG15, SSSP, PSlibrary, and GBRV all distribute in UPF. Parsing complexity is moderate — DFTK.jl's implementation is ~269 lines of Julia, a good reference. Use Rust's `quick-xml` or a custom streaming parser rather than a strict XML library.

**PSP8** is ABINIT's native norm-conserving format, designed by Hamann for ONCVPSP output. It uses plain-text fixed-format data on linear radial grids (not logarithmic), which provides better accuracy at large cutoffs by avoiding aliasing noise. The format is trivial to parse and is PseudoDojo's native output format. **PAW-XML** (version 0.7) is a true, schema-compliant XML format used by GPAW, ABINIT, and AtomPAW. It cleanly stores all PAW quantities (AE/pseudo partial waves, projectors, core densities, shape functions) in Hartree atomic units. This is the natural format for PAW implementation.

**PSML** (García et al., Comput. Phys. Commun. 227, 51, 2018) is the cleanest XML format with a formal RELAX-NG schema, supporting NC pseudopotentials for SIESTA and ABINIT. **POTCAR** (VASP) is proprietary and cannot be redistributed or officially documented. Despite this, VASP's PAW potentials are the most widely referenced in the literature and participate in all major benchmarks.

**For a new Rust code**: implement UPF2 first (largest ecosystem), PSP8 second (simplest parsing, PseudoDojo native), and PAW-XML when adding PAW support. Rust's `quick-xml` or `roxmltree` crates handle the XML formats well.

---

## Implementation architecture: from Kleinman-Bylander to GPU kernels

The **Kleinman-Bylander transformation** (Phys. Rev. Lett. 48, 1425, 1982) converts semi-local pseudopotentials into separable, fully nonlocal form: V_NL|ψ⟩ = Σ_{l,m} |χ_{lm}⟩ E_l ⟨χ_{lm}|ψ⟩. This factorizes into inner products (projections p_{lm} = ⟨χ_{lm}|ψ⟩) followed by accumulation, reducing cost from O(N²_PW) to O(N_PW × N_proj × N_bands). In reciprocal space, projections involve structure factors and can be cast as **three matrix-matrix multiplications (ZGEMM)**, directly amenable to BLAS/cuBLAS acceleration. For ONCV with multiple projectors per channel, the formalism generalizes with a D_mn matrix and multiple β projectors.

**Reciprocal-space projector application** is standard and preferred for small-to-medium systems with many k-points: projectors χ_l(|k+G|) are precomputed on the G-vector grid, and application costs O(N_PW × N_proj × N_bands) per k-point. **Real-space projection** (King-Smith, Payne, Lin, Phys. Rev. B 44, 13063, 1991) confines projectors to spheres of radius ~1.5–2.0 × r_c around each atom, reducing scaling from O(N³) to O(N²) for large systems. The crossover occurs at roughly 50–100 atoms. For a new code, implement reciprocal-space first; add real-space optimization later for large-system capability.

**PAW implementation** adds five major components beyond NCPP. Augmentation charges Q_αβ(r) = φ*_α φ_β − φ̃*_α φ̃_β represent the one-center charge difference, requiring evaluation on radial grids within PAW spheres. Compensation charges n̂(r) restore correct multipole moments on the plane-wave grid using shape functions (Gauss, sinc, or Bessel); their spheres must not overlap. VASP uses a **double-grid technique** with compensation charges on a finer FFT grid controlled by ENAUG. The total energy decomposes into three terms: Ẽ (plane-wave grid, soft quantities + compensation), E¹ (AE one-center on radial grids), and −Ẽ¹ (PS one-center on radial grids). The overlap operator S = 1 + Σ |p̃_α⟩ q_αβ ⟨p̃_β| transforms the Kohn-Sham equations into a generalized eigenvalue problem requiring S-orthogonal eigensolvers.

**GPU acceleration** targets three operations. 3D FFTs use cuFFT, potentially batching multiple bands to saturate GPU throughput. Nonlocal pseudopotential application maps to batched ZGEMM via cuBLAS — the projection, D_ij application, and accumulation steps are all matrix-matrix multiplies. QE (since v6.5) accelerates all major PP operations via OpenACC + CUDA Fortran, achieving 2–20× speedups. A recent approach by Dal Corso et al. (arXiv:2412.01695, 2024) loads **multiple k-points simultaneously** on GPU for metals with small unit cells, achieving significant additional speedups. GPU memory (16–80 GB per card) is the limiting factor: wavefunctions for all bands at one k-point consume N_PW × N_bands × 16 bytes.

**Ghost states** — spurious eigenstates of the Kleinman-Bylander operator with incorrect nodal structure — are the primary numerical stability concern. Detection methods include Bessel function diagonalization (atom-in-a-box), logarithmic derivative analysis (ghost states appear as peaks only in pseudo, not all-electron, derivatives), and SCF convergence monitoring. The ONCV two-projector approach is highly effective at suppressing ghosts compared to single-projector KB. Other remedies include changing the local channel and adjusting cutoff radii. PseudoDojo and SSSP both include automated ghost-state detection in their validation protocols.

---

## Emerging directions: from ML-optimized potentials to functional consistency

Several cutting-edge developments are reshaping the pseudopotential landscape. **Algorithmic differentiation for pseudopotential fitting** (Herbst, Levitt et al., npj Comput. Mater., 2025) uses the DFTK code's automatic differentiation capability to optimize pseudopotential parameters directly against bulk DFT observables via gradient-based methods. Currently demonstrated only for lithium, this framework could eventually enable systematically optimal pseudopotentials trained against arbitrary target properties.

**The meta-GGA consistency problem** is now largely resolved. SCAN's isoorbital indicator causes XC potential divergences that prevent reliable PP generation (Bartók and Yates, J. Chem. Phys. 150, 161101, 2019). The **r2SCAN** functional (Furness et al., J. Phys. Chem. Lett. 11, 8208, 2020) restores SCAN's accuracy while enabling smooth PP construction. In practice, **most SCAN/r2SCAN calculations still use PBE-generated pseudopotentials**, which Borlido et al. (J. Chem. Theory Comput. 16, 3620, 2020) showed introduces average errors of ~0.1 eV in band gaps — acceptable for most applications but problematic for d-electron systems where errors can exceed 1 eV. Jürg Hutter generated full GTH pseudopotential sets for SCAN and PBE0 for CP2K (github.com/juerghutter/GTH), currently the most complete functionally-consistent meta-GGA PP tables. PseudoDojo has announced but not yet released r2SCAN tables.

A provocative 2025 preprint by Ji, Lin, Ren, and He (arXiv:2505.07269) proposes **"crossed pseudopotential" calculations** — using hybrid-functional-generated PPs with PBE calculations. For 54 monovalent-Cu compounds, conventional PBE wrongly predicted ~25% as metals; crossed PPs reduced mean relative band gap error from 80% to 20%, outperforming even HSE06. This requires community validation but suggests pseudopotential quality may matter more than functional choice for certain systems.

**Correlation-consistent effective core potentials (ccECPs)** from the Mitas group bridge quantum chemistry and plane-wave DFT. The ccECP-soft variants (Kincaid et al., J. Chem. Phys. 157, 174307, 2022) are specifically adapted for plane-wave calculations with cutoffs ≤400 Ry for 3d transition metals. Coverage is expanding to lanthanides and heavy elements (Zhou et al., J. Chem. Phys. 160, 084302, 2024; Madany et al., J. Chem. Phys. 163, 114108, 2025), hosted at pseudopotentiallibrary.org. These potentials enable consistent calculations across Gaussian-basis molecular and plane-wave periodic codes — increasingly important for multi-scale workflows.

---

## How the major codes differ and what users expect

**VASP** uses PAW exclusively for production work, with proprietary POTCAR files providing one recommended potential per element (plus _s, _h, _pv, _sv, _GW variants). Its user experience is the benchmark: automatic cutoff suggestions via ENMAX, clear documentation of semicore choices, and the Materials Project's standardized POTCAR recommendations. VASP achieves ~0.3 meV/atom in Delta-test benchmarks. SOC uses on-site one-center corrections within PAW spheres.

**Quantum ESPRESSO** supports NC, USPP, and PAW via UPF format, but forces users to choose among many PP options — a frequent pain point. The SSSP library partially addresses this by curating the best PP per element. QE's open-source nature means all PP libraries are freely redistributable. GPU acceleration (v6.5+) covers all major PP operations.

**CASTEP** uniquely offers **on-the-fly pseudopotential generation (OTFG)** from compact specification strings, achieving 0.4 meV/atom Delta-test accuracy for ultrasoft and 1.1 meV/atom for norm-conserving. Users never deal with PP files directly — arguably the best user experience. This approach also ensures automatic functional consistency.

**ABINIT** is tightly integrated with PseudoDojo and JTH PAW tables, supporting psp8 and PAW-XML formats. Its **LibPAW** library (Fortran, extracted from ABINIT sources) provides reusable PAW data structures and operations, though for a Rust code it is better used as algorithmic reference than as a direct dependency via FFI.

---

## Concrete implementation roadmap for the Rust code

**Phase 1 (MVP, ~2–3 months): ONCV with UPF2.** Parse UPF2 files for NC pseudopotentials. Extract radial grid, local potential, beta projectors, D_ij matrix, NLCC density, pseudo-atomic wavefunctions. Interpolate radial functions to reciprocal space via Fourier-Bessel transform. Apply Kleinman-Bylander nonlocal PP. Target PseudoDojo NC (standard) as the primary library. Add PSP8 parsing for direct PseudoDojo access. Ship with element-to-PP mapping metadata for zero-configuration defaults. DFTK.jl's PspUpf.jl (~269 lines of Julia) and PseudoPotentialData.jl are excellent references.

**Phase 2 (~2–3 months): Expand NC + prepare PAW.** Add SG15 compatibility, automatic cutoff recommendations from PP metadata, spin-polarized calculations, and begin PAW-XML parsing. Design data structures to accommodate augmentation charges and one-center corrections from the start.

**Phase 3 (~3–6 months): PAW + SOC.** Implement augmentation charges, compensation charges, PAW overlap operator, three-term energy decomposition, and modified eigensolvers. Add spin-orbit coupling with fully-relativistic pseudopotentials (j-dependent projectors for NC; on-site corrections for PAW). Target JTH PAW table compatibility.

**Phase 4: Advanced capabilities.** DFT+U with occupation matrix control (essential for actinides), hybrid functional support, meta-GGA with kinetic energy density, DMFT projector generation, and real-space projector optimization for large systems.

Three format priorities emerge clearly. **UPF2 is non-negotiable** — it accesses the QE ecosystem (PseudoDojo, SSSP, SG15, GBRV, PSlibrary). **PSP8** is the simplest format to parse and is PseudoDojo's native output. **PAW-XML** is the cleanest, most portable PAW format. PSML is worth supporting eventually for SIESTA interoperability. POTCAR cannot be supported due to licensing restrictions, but users migrating from VASP should find equivalent accuracy through PseudoDojo PAW datasets.

## Conclusion

The pseudopotential ecosystem has reached remarkable maturity. The 2024 Bosoni et al. verification study demonstrates that modern codes and libraries converge to within ~1 meV/atom of all-electron results across 960 equations of state spanning the entire periodic table. For a new Rust-based code, the strategic path is clear: **ONCV via UPF2/PSP8 for rapid capability, PAW via PAW-XML for competitive accuracy, PseudoDojo as the primary library.** The 2024 actinide extension makes even UO₂ calculations accessible without proprietary potentials. Three emerging trends deserve architectural attention: functional-consistent pseudopotentials (design PP data structures to store functional metadata), ML-optimized generation (the AD-DFPT framework suggests future PP fitting will be gradient-based), and the convergence of quantum chemistry ECPs with plane-wave potentials (ccECP-soft variants already bridge this gap). The code that best automates PP selection — approaching CASTEP's OTFG simplicity while offering VASP's accuracy — will win users from both camps.