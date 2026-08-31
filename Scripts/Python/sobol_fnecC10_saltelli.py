#!/usr/bin/env python3
"""
Canonical Sobol GSA for the 10-input trait-resolved fnecC model.

FINAL MODEL CHOICES
-------------------
1. MS and GS are sampled independently.
2. NO GS >= MS constraint is imposed.
3. rB correction IS used here because this is the trait-resolved model.
4. rB is sampled as a MOLAR GlcN:MurA ratio and converted to a MASS
   ratio before correcting mass-based soil GlcN concentrations.
5. Trait-resolved equations are:

       CFB = cB / (MGP*fGP + MGN*(1-fGP))
       CFF = cF / GF

       NB = CFB * MS

       rB_mass = rB_molar * MW_GLCN / MW_MURA
       GS_corr = max(GS - rB_mass*MS, 0)
       NF = CFF * GS_corr

       fnecC_raw = 100 * (NB + NF) / S

6. Sobol indices are calculated for both the raw and physically
   bounded responses:

       fnecC_bounded = min(fnecC_raw, 100)

   The raw response is the primary main-text sensitivity target.
   The bounded response is retained as a robustness analysis and is
   reported alongside the raw analysis in Appendix C.

7. The complete Sobol sampling matrix is retained; no Sobol rows are
   rejected or filtered.

INPUT DISTRIBUTIONS
-------------------
Empirical marginal distributions are reconstructed from the canonical
LHS output:

    data/derived/lhs/fnecC10_lhs.csv

The Sobol analysis varies all 10 physical inputs independently through
their empirical marginals. The row-wise relationship between MS, GS,
and S in the LHS source file is NOT imposed on the Sobol matrix.

RUN FROM REPOSITORY ROOT
------------------------
python -m analysis.Python.sobol_fnecC10_saltelli

OUTPUTS
-------
data/derived/gsa/sobol_fnecC10_raw_saltelli_summary.csv
data/derived/gsa/sobol_fnecC10_bounded_saltelli_summary.csv
data/derived/gsa/sobol_fnecC10_raw_vs_bounded.csv
data/derived/gsa/sobol_fnecC10_diagnostics.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from SALib.sample.sobol import sample as sobol_sample
from SALib.analyze.sobol import analyze as sobol_analyze

from . import model_definitions as md


# =====================================================================
# Helpers
# =====================================================================

def is_power_of_two(n):
    """Return True when n is a positive power of two."""
    return n > 0 and (n & (n - 1) == 0)


def make_invq(values):
    """
    Construct an empirical inverse-CDF mapping from finite values.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    x = np.sort(x)

    if x.size < 50:
        raise ValueError(
            f"Too few finite values to build empirical inverse CDF: n={x.size}"
        )

    u_grid = np.linspace(
        0.0,
        1.0,
        x.size,
        endpoint=True,
    )

    def invq(u):
        u = np.asarray(u, dtype=float)
        u = np.clip(u, 1e-12, 1.0 - 1e-12)
        return np.interp(u, u_grid, x)

    return invq


# =====================================================================
# Main analysis
# =====================================================================

def run_sobol(
    n_base=8192,
    outdir=None,
    lhs_path=None,
    seed=1234,
    num_resamples=100,
    save_draws=True,
):
    """
    Run the canonical 10-input trait-resolved fnecC Sobol analysis.
    """

    if not is_power_of_two(n_base):
        raise ValueError(
            f"n_base={n_base} is not a power of two. "
            "Use 1024, 2048, 4096, 8192, 16384, etc."
        )

    if outdir is None:
        outdir = (
            md.REPO_ROOT
            / "data"
            / "derived"
            / "gsa"
        )
    else:
        outdir = Path(outdir)

    if lhs_path is None:
        lhs_path = (
            md.REPO_ROOT
            / "data"
            / "derived"
            / "lhs"
            / "fnecC10_lhs.csv"
        )
    else:
        lhs_path = Path(lhs_path)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not lhs_path.exists():
        raise FileNotFoundError(
            f"Canonical trait-resolved LHS file not found: {lhs_path}"
        )

    print("=" * 72)
    print("[fnecC10] Trait-resolved Sobol global sensitivity analysis")
    print("=" * 72)
    print("[fnecC10] LHS source        :", lhs_path)
    print("[fnecC10] output directory :", outdir)
    print("[fnecC10] base N           :", n_base)
    print("[fnecC10] seed             :", seed)
    print("[fnecC10] GS >= MS         : NO")
    print("[fnecC10] rB correction    : YES (trait-resolved only)")
    print("[fnecC10] rB units         : molar input -> mass correction")
    print("[fnecC10] Sobol targets    : raw and bounded at 100%")
    print("[fnecC10] Main-text target : raw, unbounded fnecC")
    print("[fnecC10] Robustness target: bounded fnecC")
    print()

    # -----------------------------------------------------------------
    # Load canonical LHS output
    # -----------------------------------------------------------------

    lhs = pd.read_csv(lhs_path)

    source_names = [
        "MS",
        "GS",
        "S",
        "cB",
        "MGP",
        "MGN",
        "fGP",
        "cF",
        "GF",
        "rB_molar",
    ]

    missing = [
        name
        for name in source_names
        if name not in lhs.columns
    ]

    if missing:
        raise ValueError(
            f"{lhs_path} missing required columns: {missing}"
        )

    invq = {
        name: make_invq(
            lhs[name].to_numpy()
        )
        for name in source_names
    }

    # -----------------------------------------------------------------
    # Independent Sobol coordinates
    # -----------------------------------------------------------------

    unit_names = [
        "u_MS",
        "u_GS",
        "u_S",
        "u_cB",
        "u_MGP",
        "u_MGN",
        "u_fGP",
        "u_cF",
        "u_GF",
        "u_rB",
    ]

    real_names = [
        "MS",
        "GS",
        "S",
        "cB",
        "MGP",
        "MGN",
        "fGP",
        "cF",
        "GF",
        "rB",
    ]

    problem = {
        "num_vars": 10,
        "names": unit_names,
        "bounds": [[0.0, 1.0]] * 10,
    }

    U = sobol_sample(
        problem,
        N=n_base,
        calc_second_order=False,
        scramble=True,
        seed=seed,
    )

    expected_rows = n_base * (10 + 2)

    if U.shape != (expected_rows, 10):
        raise RuntimeError(
            "Unexpected Sobol sample shape. "
            f"Expected {(expected_rows, 10)}, received {U.shape}."
        )

    print(
        "[fnecC10] model evaluations:",
        len(U)
    )

    # -----------------------------------------------------------------
    # Map each coordinate independently to its empirical marginal.
    #
    # IMPORTANT:
    # There is intentionally NO GS >= MS transformation here.
    # -----------------------------------------------------------------

    params = pd.DataFrame({
        name: invq[name](U[:, i])
        for i, name in enumerate(source_names)
    })

    MS = params["MS"].to_numpy(dtype=float)
    GS = params["GS"].to_numpy(dtype=float)
    S = params["S"].to_numpy(dtype=float)

    cB = params["cB"].to_numpy(dtype=float)
    MGP = params["MGP"].to_numpy(dtype=float)
    MGN = params["MGN"].to_numpy(dtype=float)
    fGP = params["fGP"].to_numpy(dtype=float)

    cF = params["cF"].to_numpy(dtype=float)
    GF = params["GF"].to_numpy(dtype=float)
    rB_molar = params["rB_molar"].to_numpy(dtype=float)

    if np.any(~np.isfinite(S)) or np.any(S <= 0):
        raise ValueError(
            "All sampled SOC values must be finite and > 0."
        )

    if np.any(~np.isfinite(GF)) or np.any(GF <= 0):
        raise ValueError(
            "All sampled fungal GlcN contents must be finite and > 0."
        )

    # -----------------------------------------------------------------
    # Trait-resolved model
    # -----------------------------------------------------------------

    CFB = md.CFB_fun(
        cB,
        MGP,
        MGN,
        fGP,
    )

    CFF = md.CFF_fun(
        cF,
        GF,
    )

    NB = md.NB_fun(
        CFB,
        MS,
    )

    # rB is a MOLAR GlcN:MurA ratio, while MS and GS are mass
    # concentrations. Convert before subtracting bacterial GlcN.
    rB_mass = md.rb_molar_to_mass(
        rB_molar
    )

    GS_corr = np.maximum(
        GS - rB_mass * MS,
        0.0,
    )

    NF = (
        np.asarray(CFF, dtype=float)
        * GS_corr
    )

    total_nec_raw = (
        np.asarray(NB, dtype=float)
        + NF
    )

    # Use the shared canonical model definition for the raw response.
    fnecC_raw = md.fnecC_trait_raw(
        MS=MS,
        GS=GS,
        S=S,
        cB=cB,
        MGP=MGP,
        MGN=MGN,
        fGP=fGP,
        cF=cF,
        GF=GF,
        rB_molar=rB_molar,
    )

    # Robustness response with explicit physical upper bound.
    fnecC_bounded = np.minimum(
        fnecC_raw,
        100.0,
    )

    if (
        np.any(~np.isfinite(fnecC_raw))
        or np.any(~np.isfinite(fnecC_bounded))
    ):
        raise RuntimeError(
            "Non-finite Sobol outputs encountered. "
            "No rows were removed."
        )

    # -----------------------------------------------------------------
    # Sobol analyses: raw and bounded responses
    #
    # Both analyses use the identical structured Sobol matrix.
    # The only difference is the output transformation.
    # -----------------------------------------------------------------

    Si_raw = sobol_analyze(
        problem,
        fnecC_raw,
        calc_second_order=False,
        num_resamples=num_resamples,
        conf_level=0.95,
        print_to_console=False,
        seed=seed,
    )

    Si_bounded = sobol_analyze(
        problem,
        fnecC_bounded,
        calc_second_order=False,
        num_resamples=num_resamples,
        conf_level=0.95,
        print_to_console=False,
        seed=seed,
    )

    summary_raw = pd.DataFrame({
        "param": real_names,
        "S1": Si_raw["S1"],
        "S1_conf": Si_raw["S1_conf"],
        "ST": Si_raw["ST"],
        "ST_conf": Si_raw["ST_conf"],
    })

    summary_bounded = pd.DataFrame({
        "param": real_names,
        "S1": Si_bounded["S1"],
        "S1_conf": Si_bounded["S1_conf"],
        "ST": Si_bounded["ST"],
        "ST_conf": Si_bounded["ST_conf"],
    })

    # -----------------------------------------------------------------
    # Save raw primary analysis
    # -----------------------------------------------------------------

    raw_summary_path = (
        outdir
        / "sobol_fnecC10_raw_saltelli_summary.csv"
    )

    summary_raw.to_csv(
        raw_summary_path,
        index=False,
    )

    np.savez(
        outdir
        / "sobol_fnecC10_raw_saltelli_Si.npz",
        **Si_raw,
    )

    # -----------------------------------------------------------------
    # Save bounded robustness analysis
    # -----------------------------------------------------------------

    bounded_summary_path = (
        outdir
        / "sobol_fnecC10_bounded_saltelli_summary.csv"
    )

    summary_bounded.to_csv(
        bounded_summary_path,
        index=False,
    )

    np.savez(
        outdir
        / "sobol_fnecC10_bounded_saltelli_Si.npz",
        **Si_bounded,
    )

    # -----------------------------------------------------------------
    # Raw-versus-bounded comparison
    # -----------------------------------------------------------------

    comparison = pd.DataFrame({
        "formulation": "trait_resolved",
        "param": real_names,

        "S1_raw": summary_raw["S1"],
        "S1_raw_conf": summary_raw["S1_conf"],

        "S1_bounded": summary_bounded["S1"],
        "S1_bounded_conf": summary_bounded["S1_conf"],

        "S1_change_bounded_minus_raw": (
            summary_bounded["S1"]
            - summary_raw["S1"]
        ),

        "ST_raw": summary_raw["ST"],
        "ST_raw_conf": summary_raw["ST_conf"],

        "ST_bounded": summary_bounded["ST"],
        "ST_bounded_conf": summary_bounded["ST_conf"],

        "ST_change_bounded_minus_raw": (
            summary_bounded["ST"]
            - summary_raw["ST"]
        ),

        "S1_rank_raw": summary_raw["S1"].rank(
            ascending=False,
            method="min",
        ),

        "S1_rank_bounded": summary_bounded["S1"].rank(
            ascending=False,
            method="min",
        ),

        "ST_rank_raw": summary_raw["ST"].rank(
            ascending=False,
            method="min",
        ),

        "ST_rank_bounded": summary_bounded["ST"].rank(
            ascending=False,
            method="min",
        ),
    })

    comparison_path = (
        outdir
        / "sobol_fnecC10_raw_vs_bounded.csv"
    )

    comparison.to_csv(
        comparison_path,
        index=False,
    )

    # -----------------------------------------------------------------
    # Diagnostics
    # -----------------------------------------------------------------

    molecular_C_feasible = md.molecular_carbon_feasible(
        MS,
        GS,
        S,
    )

    S1_rank_spearman = summary_raw["S1"].rank().corr(
        summary_bounded["S1"].rank()
    )

    ST_rank_spearman = summary_raw["ST"].rank().corr(
        summary_bounded["ST"].rank()
    )

    diagnostics = pd.DataFrame({
        "metric": [
            "n_base",
            "n_model_evaluations",
            "seed",
            "prop_GS_lt_MS",
            "prop_molecular_C_infeasible",
            "n_molecular_C_infeasible",
            "raw_median",
            "raw_p99",
            "raw_max",
            "prop_raw_gt_100",
            "n_raw_gt_100",
            "raw_variance",
            "bounded_median",
            "bounded_max",
            "prop_at_upper_bound",
            "bounded_variance",
            "variance_ratio_bounded_to_raw",
            "S1_rank_spearman_raw_vs_bounded",
            "ST_rank_spearman_raw_vs_bounded",
            "rB_molar_min",
            "rB_molar_max",
            "rB_mass_min",
            "rB_mass_max",
        ],
        "value": [
            n_base,
            len(fnecC_raw),
            seed,
            np.mean(GS < MS),
            np.mean(~molecular_C_feasible),
            np.sum(~molecular_C_feasible),
            np.median(fnecC_raw),
            np.quantile(fnecC_raw, 0.99),
            np.max(fnecC_raw),
            np.mean(fnecC_raw > 100.0),
            np.sum(fnecC_raw > 100.0),
            np.var(fnecC_raw, ddof=1),
            np.median(fnecC_bounded),
            np.max(fnecC_bounded),
            np.mean(
                np.isclose(
                    fnecC_bounded,
                    100.0,
                )
            ),
            np.var(fnecC_bounded, ddof=1),
            (
                np.var(fnecC_bounded, ddof=1)
                / np.var(fnecC_raw, ddof=1)
            ),
            S1_rank_spearman,
            ST_rank_spearman,
            np.min(rB_molar),
            np.max(rB_molar),
            np.min(rB_mass),
            np.max(rB_mass),
        ],
    })

    diagnostics.to_csv(
        outdir
        / "sobol_fnecC10_diagnostics.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # Save exact evaluated draws
    # -----------------------------------------------------------------

    if save_draws:

        pd.DataFrame(
            U,
            columns=unit_names,
        ).to_csv(
            outdir
            / "sobol_fnecC10_draws_uniform.csv",
            index=False,
        )

        params_out = params.copy()

        params_out["rB_mass"] = rB_mass
        params_out["GS_corr"] = GS_corr
        params_out["CFB"] = CFB
        params_out["CFF"] = CFF

        params_out.to_csv(
            outdir
            / "sobol_fnecC10_draws_params.csv",
            index=False,
        )

        pd.DataFrame({
            "NB": NB,
            "NF": NF,
            "totalNec_raw": total_nec_raw,
            "fnecC_raw": fnecC_raw,
            "fnecC_bounded": fnecC_bounded,
        }).to_csv(
            outdir
            / "sobol_fnecC10_outputs.csv",
            index=False,
        )

    # -----------------------------------------------------------------
    # Metadata
    # -----------------------------------------------------------------

    metadata = pd.DataFrame({
        "item": [
            "analysis",
            "n_base",
            "num_vars",
            "sampler",
            "seed",
            "sampling_space",
            "marginal_mapping",
            "GS_MS_constraint",
            "rB_correction",
            "rB_units",
            "rB_conversion",
            "NB_equation",
            "NF_equation",
            "raw_output",
            "sobol_targets",
            "main_text_target",
            "robustness_target",
            "appendix_location",
            "row_filtering",
        ],
        "details": [
            "Canonical 10-input trait-resolved fnecC Sobol GSA",
            str(n_base),
            "10",
            "SALib.sample.sobol.sample",
            str(seed),
            "Independent unit-hypercube Sobol coordinates",
            (
                "Each coordinate transformed independently through "
                "its empirical marginal from fnecC10_lhs.csv"
            ),
            "None; GS is NOT constrained to be >= MS",
            "Applied only in trait-resolved model",
            "rB sampled as molar GlcN:MurA ratio",
            (
                "rB_mass = rB_molar * MW_GLCN / MW_MURA"
            ),
            "NB = CFB * MS",
            (
                "NF = CFF * max(GS - rB_mass*MS, 0)"
            ),
            "fnecC_raw = 100*(NB + NF)/S",
            (
                "Raw fnecC and "
                "fnecC_bounded = min(fnecC_raw, 100)"
            ),
            "Raw, unbounded fnecC",
            "fnecC bounded at 100% of SOC",
            "Raw-versus-bounded comparison reported in Appendix C",
            "None; complete Sobol matrix retained",
        ],
    })

    metadata.to_csv(
        outdir
        / "sobol_fnecC10_metadata.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # Console summary
    # -----------------------------------------------------------------

    print()
    print("[fnecC10] Raw Sobol indices (primary main-text target)")
    print(
        summary_raw.sort_values(
            "ST",
            ascending=False,
        ).to_string(
            index=False
        )
    )

    print()
    print("[fnecC10] Bounded Sobol indices (robustness analysis)")
    print(
        summary_bounded.sort_values(
            "ST",
            ascending=False,
        ).to_string(
            index=False
        )
    )

    print()
    print(
        "[fnecC10] GS < MS diagnostic:",
        f"{100 * np.mean(GS < MS):.2f}%"
    )

    print(
        "[fnecC10] molecular-C input infeasible:",
        f"{100 * np.mean(~molecular_C_feasible):.4f}%"
    )

    print(
        "[fnecC10] raw fnecC > 100%:",
        f"{100 * np.mean(fnecC_raw > 100.0):.2f}%"
    )

    print(
        "[fnecC10] ST rank Spearman, raw vs bounded:",
        f"{ST_rank_spearman:.4f}"
    )

    print(
        "[fnecC10] saved:",
        raw_summary_path
    )

    print(
        "[fnecC10] saved:",
        bounded_summary_path
    )

    print(
        "[fnecC10] saved:",
        comparison_path
    )

    # Primary return object now matches the main-text sensitivity target.
    return summary_raw, Si_raw, diagnostics


# =====================================================================
# Command-line interface
# =====================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Canonical Sobol GSA for the 10-input trait-resolved fnecC model."
        )
    )

    parser.add_argument(
        "--n_base",
        type=int,
        default=8192,
    )

    parser.add_argument(
        "--outdir",
        default=None,
    )

    parser.add_argument(
        "--lhs_path",
        default=None,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
    )

    parser.add_argument(
        "--num_resamples",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--no_save_draws",
        action="store_true",
    )

    args = parser.parse_args()

    run_sobol(
        n_base=args.n_base,
        outdir=args.outdir,
        lhs_path=args.lhs_path,
        seed=args.seed,
        num_resamples=args.num_resamples,
        save_draws=(
            not args.no_save_draws
        ),
    )