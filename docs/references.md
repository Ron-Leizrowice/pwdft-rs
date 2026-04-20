# Key References

Canonical papers and textbooks cited in this codebase.

## Textbooks

- **Martin, R. M.** — *Electronic Structure: Basic Theory and Practical
  Methods*, 2nd ed., Cambridge University Press (2020).
  Chapters 11-13 cover plane-wave pseudopotential DFT.

- **Ashcroft & Mermin** — *Solid State Physics*, Saunders College (1976).
  Chapter 17: Thomas-Fermi screening (used in Kerker preconditioning).

## Review Papers

- **Payne, Teter, Allan, Arias, Joannopoulos** — *Iterative minimization
  techniques for ab initio total-energy calculations*, Rev. Mod. Phys.
  **64**, 1045 (1992). The canonical review of plane-wave PP DFT.

## Exchange-Correlation

- **Slater** — Phys. Rev. **81**, 385 (1951). Slater exchange.

- **Ceperley & Alder** — Phys. Rev. Lett. **45**, 566 (1980).
  Quantum Monte Carlo data for the homogeneous electron gas.

- **Perdew & Zunger** — Phys. Rev. B **23**, 5048 (1981).
  Parametrization of Ceperley-Alder data. Table I has all our LDA parameters.

## Pseudopotentials

- **Kleinman & Bylander** — Phys. Rev. Lett. **48**, 1425 (1982).
  Separable form for non-local pseudopotentials (KB projectors).

- **Louie, Froyen, Cohen** — Phys. Rev. B **26**, 1738 (1982).
  Nonlinear core correction (NLCC).

## Smearing

- **Methfessel & Paxton** — Phys. Rev. B **40**, 3616 (1989).
  Hermite polynomial smearing.

- **Marzari, Vanderbilt, De Vita, Payne** — Phys. Rev. Lett. **82**, 3296 (1999).
  Cold smearing (Marzari-Vanderbilt).

## Methodology

- **Kresse & Furthmuller** — Phys. Rev. B **54**, 11169 (1996).
  VASP methodology (Davidson eigensolver, Pulay mixing).

## Online Resources

- **ABINIT theory docs** — <https://docs.abinit.org/theory/pseudopotentials/>
- **Theoretical Physics Reference** — <https://www.theoretical-physics.com/dev/quantum/dft.html>
- **MAVENs DFT Notes** — <https://mavens-group.github.io/dft-notes/08-SCF.html>
