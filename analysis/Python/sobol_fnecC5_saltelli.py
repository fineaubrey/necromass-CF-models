#!/usr/bin/env python3
"""
Canonical Sobol GSA for the five-input composite fnecC model.

FINAL MODEL CHOICES
-------------------
1. MS and GS are sampled independently.
2. NO GS >= MS constraint is imposed.
3. NO rB correction is used in the composite formulation.
4. Composite equations are:

       NB = CFB * MS
       NF = CFF * GS
       fnecC_raw = 100 * (NB + NF) / S

5. Sobol indices are calculated for both the raw and physically
   bounded responses:

       fnecC_bounded = min(fnecC_raw, 100)

   The bounded response remains the canonical main-text target so
   existing Figure 3 workflows remain unchanged. Raw indices are
   retained for comparison in Appendix D.

6. The complete Sobol sampling matrix is retained; no Sobol rows are
   rejected or filtered.

INPUT DISTRIBUTIONS
-------------------
Empirical marginal distributions are reconstructed from the canonical
LHS output:

    data/derived/lhs/fnecC5_lhs.csv

The LHS source file may itself reflect feasibility filtering used during
uncertainty propagation, but the Sobol analysis varies the five physical
inputs independently through their empirical marginals.

RUN FROM REPOSITORY ROOT
------------------------
python -m analysis.Python.sobol_fnecC5_saltelli

OUTPUTS
-------
data/derived/gsa/sobol_fnecC5_saltelli_summary.csv
data/derived/gsa/sobol_fnecC5_raw_saltelli_summary.csv
data/derived/gsa/sobol_fnecC5_raw_vs_bounded.csv
data/derived/gsa/sobol_fnecC5_diagnostics.csv
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
    n_base=4096,
    outdir=None,
    lhs_path=None,
    seed=1234,
    num_resamples=100,
    save_draws=True,
):
    """
    Run the canonical five-input composite fnecC Sobol analysis.
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
            / "fnecC5_lhs.csv"
        )
    else:
        lhs_path = Path(lhs_path)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not lhs_path.exists():
        raise FileNotFoundError(
            f"Canonical composite LHS file not found: {lhs_path}"
        )

    print("=" * 72)
    print("[fnecC5] Composite Sobol global sensitivity analysis")
    print("=" * 72)
    print("[fnecC5] LHS source        :", lhs_path)
    print("[fnecC5] output directory :", outdir)
    print("[fnecC5] base N           :", n_base)
    print("[fnecC5] seed             :", seed)
    print("[fnecC5] GS >= MS         : NO")
    print("[fnecC5] rB correction    : NO")
    print("[fnecC5] Sobol targets    : raw and bounded at 100%")
    print("[fnecC5] Main-text target : bounded at 100%")
    print()

    # -----------------------------------------------------------------
    # Load canonical LHS output
    # -----------------------------------------------------------------

    lhs = pd.read_csv(lhs_path)

    parameter_names = [
        "MS",
        "GS",
        "S",
        "CFB",
        "CFF",
    ]

    missing = [
        name
        for name in parameter_names
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
        for name in parameter_names
    }

    # -----------------------------------------------------------------
    # Independent Sobol coordinates
    # -----------------------------------------------------------------

    unit_names = [
        "u_MS",
        "u_GS",
        "u_S",
        "u_CFB",
        "u_CFF",
    ]

    problem = {
        "num_vars": 5,
        "names": unit_names,
        "bounds": [[0.0, 1.0]] * 5,
    }

    U = sobol_sample(
        problem,
        N=n_base,
        calc_second_order=False,
        scramble=True,
        seed=seed,
    )

    expected_rows = n_base * (5 + 2)

    if U.shape != (expected_rows, 5):
        raise RuntimeError(
            "Unexpected Sobol sample shape. "
            f"Expected {(expected_rows, 5)}, received {U.shape}."
        )

    print(
        "[fnecC5] model evaluations:",
        len(U)
    )

    # -----------------------------------------------------------------
    # Map each coordinate independently to its empirical marginal
    #
    # IMPORTANT:
    # There is intentionally NO GS >= MS transformation here.
    # -----------------------------------------------------------------

    params = pd.DataFrame({
        name: invq[name](U[:, i])
        for i, name in enumerate(parameter_names)
    })

    MS = params["MS"].to_numpy(dtype=float)
    GS = params["GS"].to_numpy(dtype=float)
    S = params["S"].to_numpy(dtype=float)
    CFB = params["CFB"].to_numpy(dtype=float)
    CFF = params["CFF"].to_numpy(dtype=float)

    if np.any(~np.isfinite(S)) or np.any(S <= 0):
        raise ValueError(
            "All sampled SOC values must be finite and > 0."
        )

    # -----------------------------------------------------------------
    # Canonical composite model
    #
    # NO rB correction.
    # NO subtraction from GS.
    # -----------------------------------------------------------------

    NB = md.NB_fun(
        CFB,
        MS,
    )

    # Composite formulation: NO bacterial GlcN correction.
    NF = (
        np.asarray(CFF, dtype=float)
        * np.asarray(GS, dtype=float)
    )

    total_nec_raw = NB + NF

    # Use the canonical shared model definition so the Sobol script
    # cannot silently diverge from model_definitions.py.
    fnecC_raw = md.fnecC_composite_raw(
        MS=MS,
        GS=GS,
        S=S,
        CFB=CFB,
        CFF=CFF,
    )

    # Physical bound used only as the Sobol response.
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
    # Both analyses use the identical structured Sobol matrix. The
    # only difference is the output transformation, which isolates the
    # effect of bounding on the variance decomposition.
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
        "param": parameter_names,
        "S1": Si_raw["S1"],
        "S1_conf": Si_raw["S1_conf"],
        "ST": Si_raw["ST"],
        "ST_conf": Si_raw["ST_conf"],
    })

    summary_bounded = pd.DataFrame({
        "param": parameter_names,
        "S1": Si_bounded["S1"],
        "S1_conf": Si_bounded["S1_conf"],
        "ST": Si_bounded["ST"],
        "ST_conf": Si_bounded["ST_conf"],
    })

    # Preserve the existing canonical filename for the bounded target
    # so downstream Figure 3 code continues to work unchanged.
    summary_path = (
        outdir
        / "sobol_fnecC5_saltelli_summary.csv"
    )

    summary_bounded.to_csv(
        summary_path,
        index=False,
    )

    np.savez(
        outdir
        / "sobol_fnecC5_saltelli_Si.npz",
        **Si_bounded,
    )

    raw_summary_path = (
        outdir
        / "sobol_fnecC5_raw_saltelli_summary.csv"
    )

    summary_raw.to_csv(
        raw_summary_path,
        index=False,
    )

    np.savez(
        outdir
        / "sobol_fnecC5_raw_saltelli_Si.npz",
        **Si_raw,
    )

    comparison = pd.DataFrame({
        "formulation": "composite",
        "param": parameter_names,
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
        / "sobol_fnecC5_raw_vs_bounded.csv"
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
        ],
    })

    diagnostics.to_csv(
        outdir
        / "sobol_fnecC5_diagnostics.csv",
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
            / "sobol_fnecC5_draws_uniform.csv",
            index=False,
        )

        params.to_csv(
            outdir
            / "sobol_fnecC5_draws_params.csv",
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
            / "sobol_fnecC5_outputs.csv",
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
            "composite_NB",
            "composite_NF",
            "raw_output",
            "sobol_targets",
            "main_text_target",
            "row_filtering",
        ],
        "details": [
            "Canonical five-input composite fnecC Sobol GSA",
            str(n_base),
            "5",
            "SALib.sample.sobol.sample",
            str(seed),
            "Independent unit-hypercube Sobol coordinates",
            (
                "Each coordinate transformed independently through "
                "its empirical marginal from fnecC5_lhs.csv"
            ),
            "None; GS is NOT constrained to be >= MS",
            "None",
            "NB = CFB * MS",
            "NF = CFF * GS",
            "fnecC_raw = 100*(NB + NF)/S",
            "Raw fnecC and fnecC_bounded = min(fnecC_raw, 100)",
            "Bounded fnecC; raw indices retained for Appendix D",
            "None; complete Sobol matrix retained",
        ],
    })

    metadata.to_csv(
        outdir
        / "sobol_fnecC5_metadata.csv",
        index=False,
    )

    print()
    print("[fnecC5] Bounded Sobol indices (canonical main-text target)")
    print(
        summary_bounded.sort_values(
            "ST",
            ascending=False,
        ).to_string(
            index=False
        )
    )

    print()
    print("[fnecC5] Raw Sobol indices (Appendix D comparison)")
    print(
        summary_raw.sort_values(
            "ST",
            ascending=False,
        ).to_string(
            index=False
        )
    )

    print()
    print(
        "[fnecC5] GS < MS diagnostic:",
        f"{100 * np.mean(GS < MS):.2f}%"
    )
    print(
        "[fnecC5] molecular-C input infeasible:",
        (
            f"{100 * np.mean(~molecular_C_feasible):.4f}%"
        )
    )
    print(
        "[fnecC5] raw fnecC > 100%:",
        f"{100 * np.mean(fnecC_raw > 100.0):.2f}%"
    )
    print(
        "[fnecC5] ST rank Spearman, raw vs bounded:",
        f"{ST_rank_spearman:.4f}"
    )
    print(
        "[fnecC5] saved:",
        summary_path
    )
    print(
        "[fnecC5] saved:",
        raw_summary_path
    )
    print(
        "[fnecC5] saved:",
        comparison_path
    )

    # Preserve the historical return signature using the bounded
    # target, which remains the canonical main-text analysis.
    return summary_bounded, Si_bounded, diagnostics


# =====================================================================
# Command-line interface
# =====================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Canonical Sobol GSA for the five-input composite fnecC model."
        )
    )

    parser.add_argument(
        "--n_base",
        type=int,
        default=4096,
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
