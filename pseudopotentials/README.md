# Pseudopotential Library

Norm-conserving pseudopotentials for plane-wave DFT calculations.

## Directory Layout

```
pseudopotentials/
├── nc/                         # Norm-conserving
│   ├── lda/                    # LDA exchange-correlation (PZ parametrization)
│   │   ├── Al.upf              # PseudoDojo ONCV standard accuracy
│   │   ├── Si.upf              # PseudoDojo ONCV standard accuracy
│   │   ├── ...                 # 70 elements from PseudoDojo
│   │   └── legacy/             # PPs from initial development/validation
│   │       ├── Si_hgh.UPF      # Hartwigsen-Goedecker-Hutter (QE distribution)
│   │       ├── Fe_dalcorso.UPF # Dal Corso NC (QE distribution)
│   │       ├── C_fhi_v1.UPF    # Fritz-Haber-Institute (UPF v1 format)
│   │       ├── Ga_oncv.UPF     # ONCV (QE GaN example)
│   │       └── N_oncv.UPF      # ONCV (QE GaN example)
│   └── pbe/                    # PBE exchange-correlation
│       ├── Al.upf
│       └── ...                 # 72 elements from PseudoDojo
└── README.md                   # This file
```

## Sources and Attribution

### PseudoDojo (nc/lda/ and nc/pbe/)

ONCV norm-conserving pseudopotentials from the PseudoDojo project,
scalar-relativistic, standard accuracy.

- **Website**: <http://www.pseudo-dojo.org>
- **Format**: UPF v2
- **Type**: Norm-conserving, optimized norm-conserving Vanderbilt (ONCV)
- **Elements**: 70 (LDA), 72 (PBE) — H through Zr plus selected heavier elements

**Citation**:

> M.J. van Setten, M. Giantomassi, E. Bousquet, M.J. Verstraete, D.R. Hamann,
> X. Gonze, G.-M. Rignanese,
> "The PseudoDojo: Training and grading a 85 element optimized norm-conserving
> pseudopotential table",
> *Computer Physics Communications* **226**, 39-54 (2018).
> DOI: [10.1016/j.cpc.2018.01.012](https://doi.org/10.1016/j.cpc.2018.01.012)

> D.R. Hamann,
> "Optimized norm-conserving Vanderbilt pseudopotentials",
> *Physical Review B* **88**, 085117 (2013).
> DOI: [10.1103/PhysRevB.88.085117](https://doi.org/10.1103/PhysRevB.88.085117)

### Legacy pseudopotentials (`nc/lda/legacy/`)

Used during initial development and validation. Referenced by existing
tests. These may differ from the PseudoDojo versions of the same elements.

| File | Type | XC | Source |
|------|------|-----|--------|
| `Si_hgh.UPF` | HGH separable | LDA (PZ) | QE distribution (`Si.pz-hgh.UPF`) |
| `Fe_dalcorso.UPF` | NC | LDA (PZ) | QE distribution (`Fe.pz-n-nc.UPF`), Dal Corso |
| `C_fhi_v1.UPF` | NC | LDA (PZ) | QE distribution, Fritz-Haber-Institute (UPF v1) |
| `Ga_oncv.UPF` | ONCV | LDA | QE EPW GaN example (`Ga_ONCV_LDA-1.0.upf`) |
| `N_oncv.UPF` | ONCV | LDA | QE EPW GaN example (`N_ONCV_LDA-1.0.upf`) |

## Supported formats

pwdft-rs currently supports:

- **UPF v2** (XML-like): Full support for norm-conserving PPs with
  local potential, KB non-local projectors, D_ij matrix, PP_RHOATOM,
  and PP_NLCC (nonlinear core correction).
- **UPF v1** (plain text): Not yet supported. C.UPF is v1 format.
- **PSP8** (ABINIT): Partial support (D_ij always zero — see proposal 18).

## Notes

- **NLCC**: Many PseudoDojo PPs include nonlinear core corrections
  (`core_correction="T"`). pwdft-rs handles these by adding the core
  charge to the valence density before XC evaluation.
- **LDA vs PBE**: Use `nc/lda/` with LDA functionals and `nc/pbe/`
  with PBE. Mixing LDA PPs with PBE XC (or vice versa) produces
  incorrect results — the PP and functional must be self-consistent.
- **Cutoff recommendations**: PseudoDojo provides suggested cutoffs
  in their database. Typical values: 30-50 Ry for first-row elements,
  40-80 Ry for transition metals.
