# Pseudopotential Library

Norm-conserving pseudopotentials for plane-wave DFT calculations.

## Directory Layout

```text
pseudopotentials/
├── nc/                         # Norm-conserving
│   ├── lda/                    # LDA exchange-correlation (PZ parametrization)
│   │   ├── Al.upf              # PseudoDojo ONCV standard accuracy
│   │   ├── Si.upf              # PseudoDojo ONCV standard accuracy
│   │   └── ...                 # 70 elements from PseudoDojo
│   └── pbe/                    # PBE exchange-correlation
│       └── ...                 # 72 elements from PseudoDojo
├── uspp/                       # Ultrasoft (SSSP efficiency, for future use)
│   └── pbe/
├── paw/                        # PAW (SSSP efficiency, for future use)
│   └── pbe/
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
>
> D.R. Hamann,
> "Optimized norm-conserving Vanderbilt pseudopotentials",
> *Physical Review B* **88**, 085117 (2013).
> DOI: [10.1103/PhysRevB.88.085117](https://doi.org/10.1103/PhysRevB.88.085117)

### SSSP — Standard Solid-State Pseudopotentials (uspp/pbe/, paw/pbe/)

Curated library from the SSSP project, selecting the best-performing
PP for each element from multiple sources (PseudoDojo, GBRV, etc.).
The "efficiency" set balances accuracy and computational cost.

- **Website**: <https://www.materialscloud.org/discover/sssp>
- **Version**: 1.3.0
- **Note**: SSSP is primarily USPP and PAW — not yet supported by pwdft-rs.
  Included for future use when ultrasoft/PAW support is added.

**Citation**:

> G. Prandini, A. Marrazzo, I.E. Castelli, N. Mounet, N. Marzari,
> "Precision and efficiency in solid-state pseudopotential calculations",
> *npj Computational Materials* **4**, 72 (2018).
> DOI: [10.1038/s41524-018-0127-2](https://doi.org/10.1038/s41524-018-0127-2)

## Supported formats

pwdft-rs currently supports:

- **UPF v2** (XML-like): Full support for norm-conserving PPs with
  local potential, KB non-local projectors, D_ij matrix, PP_RHOATOM,
  and PP_NLCC (nonlinear core correction).

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
