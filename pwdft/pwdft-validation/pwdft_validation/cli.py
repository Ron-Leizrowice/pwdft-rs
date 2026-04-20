"""pwdft-validate CLI — cyclopts-based sub-command dispatcher."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import cyclopts

from pwdft_validation import density as _density
from pwdft_validation import energy as _energy
from pwdft_validation import fermi as _fermi
from pwdft_validation.diagnostics import ef_shift as _ef_shift
from pwdft_validation.diagnostics import pbe_refs as _pbe_refs
from pwdft_validation.diagnostics import residual_scan as _residual_scan
from pwdft_validation.paths import CSV_REF_DIR, PSEUDO_DIR, QE_REF_DIR
from pwdft_validation.reference import beta as _beta
from pwdft_validation.reference import dij as _dij
from pwdft_validation.reference import hamiltonian as _hamiltonian
from pwdft_validation.reference import nlcc as _nlcc
from pwdft_validation.reference import sad as _sad
from pwdft_validation.reference import vloc as _vloc

app = cyclopts.App(
    name="pwdft-validate",
    help="Validation and reference-data tools for pwdft-rs.",
)

# ---------------------------------------------------------------------------
# reference sub-app
# ---------------------------------------------------------------------------

reference = cyclopts.App(name="reference", help="Generate reference CSV pin files.")
app.command(reference)


@reference.command(name="vloc")
def reference_vloc(*, heavy: bool = False) -> None:
    """V_local(G) reference.

    Without --heavy: Si FCC first-20-shell V_local(G) → vloc_g_si_reference.csv.
    With --heavy:    V_local(G=0) for 11 VGCH-scope elements → vgch_vloc_heavy.csv.
    """
    pp_dir = PSEUDO_DIR / "nc" / "lda"
    if heavy:
        raise SystemExit(_vloc.generate_heavy(pp_dir, CSV_REF_DIR / "vgch_vloc_heavy.csv"))
    raise SystemExit(_vloc.generate_si(pp_dir, CSV_REF_DIR / "vloc_g_si_reference.csv"))


@reference.command(name="beta")
def reference_beta(*, heavy: bool = False) -> None:
    """KB projector form factors β_l(q).

    Without --heavy: Si FCC (VGCMP Phase 2) → beta_q_si_reference.csv.
    With --heavy:    11 VGCH-scope elements → vgch_beta_l_heavy.csv.
    """
    pp_dir = PSEUDO_DIR / "nc" / "lda"
    if heavy:
        raise SystemExit(_beta.generate_heavy(pp_dir, CSV_REF_DIR / "vgch_beta_l_heavy.csv"))
    raise SystemExit(_beta.generate_si(pp_dir, CSV_REF_DIR / "beta_q_si_reference.csv"))


@reference.command(name="dij")
def reference_dij() -> None:
    """D_ij KB coupling matrix for Si (VGCMP Phase 3) → dij_si_reference.csv."""
    raise SystemExit(_dij.generate(PSEUDO_DIR / "nc" / "lda", CSV_REF_DIR / "dij_si_reference.csv"))


@reference.command(name="hamiltonian")
def reference_hamiltonian() -> None:
    """Assembled H diagonal for Si at k=Γ (VGCMP Phase 4) → vgcmp_phase4_reference.csv."""
    raise SystemExit(_hamiltonian.generate(PSEUDO_DIR / "nc" / "lda", CSV_REF_DIR / "vgcmp_phase4_reference.csv"))


@reference.command(name="sad")
def reference_sad() -> None:
    """SAD initial density for 7 VGCH systems → vgch_sad_heavy.csv."""
    raise SystemExit(_sad.generate(PSEUDO_DIR / "nc" / "lda", CSV_REF_DIR / "vgch_sad_heavy.csv"))


@reference.command(name="nlcc")
def reference_nlcc() -> None:
    """ρ_core(G) NLCC core-density for Si/Fe/Cu/Mn → rho_core_g_reference.csv."""
    raise SystemExit(_nlcc.generate(PSEUDO_DIR / "nc" / "lda", CSV_REF_DIR / "rho_core_g_reference.csv"))


# ---------------------------------------------------------------------------
# energy sub-app
# ---------------------------------------------------------------------------

energy = cyclopts.App(name="energy", help="Parse and compare QE energy decompositions.")
app.command(energy)


@energy.command(name="components")
def energy_components() -> None:
    """VGC5 per-component energies for Si + Fe → vgc5_qe_*.csv."""
    raise SystemExit(_energy.extract_components(QE_REF_DIR, CSV_REF_DIR))


@energy.command(name="trace")
def energy_trace(
    out: Annotated[Path, cyclopts.Parameter(help="Output CSV path.")] = CSV_REF_DIR / "vgch2_per_term_trace.csv",
) -> None:
    """VGCH-2 8-system per-term trace from QE output files."""
    raise SystemExit(_energy.extract_trace(QE_REF_DIR, out))


@energy.command(name="join")
def energy_join(
    pwdft_csv: Annotated[Path, cyclopts.Parameter(help="pwdft per-term CSV (from Tier-2 test).")] = Path(
        "target/tmp/vgch2_per_term_trace_pwdft.csv"
    ),
) -> None:
    """Join QE + pwdft per-term CSVs and print delta ranking."""
    qe_csv = CSV_REF_DIR / "vgch2_per_term_trace.csv"
    raise SystemExit(_energy.join_traces(qe_csv, pwdft_csv))


# ---------------------------------------------------------------------------
# fermi command
# ---------------------------------------------------------------------------


@app.command(name="fermi")
def fermi_command(
    out: Annotated[Path, cyclopts.Parameter(help="Output CSV path.")] = CSV_REF_DIR / "vgch2c_fermi_reference.csv",
) -> None:
    """Fermi-level bisection reference for Cu/Fe/NaCl/C/MgO → vgch2c_fermi_reference.csv."""
    raise SystemExit(_fermi.compute(QE_REF_DIR, out))


# ---------------------------------------------------------------------------
# density command
# ---------------------------------------------------------------------------


@app.command(name="density")
def density_command(
    rho: Annotated[Path, cyclopts.Parameter(help="Path to QE charge-density.dat.")],
    out: Annotated[Path, cyclopts.Parameter(help="Output .bin path (VGCH2BIN bundle).")],
    *,
    verbose: bool = False,
) -> None:
    """Parse QE Fortran binary charge-density.dat → VGCH2BIN flat binary."""
    raise SystemExit(_density.parse(rho, out, verbose=verbose))


# ---------------------------------------------------------------------------
# diag sub-app
# ---------------------------------------------------------------------------

diag = cyclopts.App(name="diag", help="Diagnostic scripts cross-checking pwdft-rs vs QE.")
app.command(diag)


@diag.command(name="ef-shift")
def diag_ef_shift() -> None:
    """Si E_F shift diagnostic: localize the ~1.35 eV rigid offset vs QE."""
    raise SystemExit(_ef_shift.run())


@diag.command(name="residual-scan")
def diag_residual_scan(
    *,
    thresholds: Annotated[list[float] | None, cyclopts.Parameter(help="conv_threshold values to sweep.")] = None,
    no_qe: bool = False,
    keep_yaml: bool = False,
    max_iter: int = 200,
) -> None:
    """PCRS per-component energy residual scan (sweeps conv_threshold)."""
    raise SystemExit(_residual_scan.run(thresholds=thresholds, no_qe=no_qe, keep_yaml=keep_yaml, max_iter=max_iter))


# ---------------------------------------------------------------------------
# pbe command
# ---------------------------------------------------------------------------


@app.command(name="pbe")
def pbe_command() -> None:
    """Extract PBE reference data from QE output logs → TOML snippet on stdout."""
    raise SystemExit(_pbe_refs.extract())


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    app()


if __name__ == "__main__":
    main()
