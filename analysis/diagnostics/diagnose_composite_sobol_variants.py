#!/usr/bin/env python3
"""
Forensic comparison of four composite-fnecC Sobol formulations.

PURPOSE
-------
This is a DIAGNOSTIC script only. It does not overwrite the production
Sobol results in data/derived/gsa/.

It tests the two methodological switches that have differed among
analysis versions:

1. Soil GlcN treatment
   A/B: GS sampled independently from its current canonical LHS marginal.
   C/D: GS sampled conditionally so that GS >= MS, reproducing the
        earlier "qGS_ge_MS" logic.

2. Sobol response
   A/C: raw, unbounded fnecC.
   B/D: fnecC bounded to 100% for Sobol analysis.

The same Sobol unit-hypercube design, random seed, MS distribution,
S distribution, CFB distribution, CFF distribution, and composite
model equation are used across all four cases.

DEFAULT COMPOSITE EQUATION
--------------------------
    NB = CFB * MS
    NF = CFF * GS
    fnecC_raw = 100 * (NB + NF) / S

This matches the current canonical composite implementation.

OPTIONAL LEGACY CORRECTION
--------------------------
To additionally test the old fixed bacterial GlcN subtraction:

    NF = CFF * max(GS - rB_mass * MS, 0)

run with:

    --legacy_rb_mass 1.16

This changes the model equation, so it should be considered a SECOND
forensic experiment rather than part of the four-case comparison.

IMPORTANT SOBOL NOTE
--------------------
In the GS>=MS cases, the physical GS value depends on both the MS and
GS Sobol coordinates. Therefore those indices are diagnostic of the
older parameterization and should not be interpreted as a classical
Sobol decomposition of independent physical MS and GS inputs.

RUN FROM REPOSITORY ROOT
------------------------
python -m analysis.Python.diagnose_composite_sobol_variants

Optional:
python -m analysis.Python.diagnose_composite_sobol_variants \
    --n_base 4096 \
    --legacy_rb_mass 1.16

OUTPUTS
-------
data/derived/gsa_diagnostics/composite_four_case_sobol.csv
data/derived/gsa_diagnostics/composite_four_case_diagnostics.csv
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
    return n > 0 and (n & (n - 1) == 0)


def make_invq(x):
    """Linear empirical inverse CDF."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = np.sort(x)

    if x.size < 50:
        raise ValueError(
            f"Too few finite values to build inverse CDF: n={x.size}"
        )

    u_grid = np.linspace(
        0.0,
        1.0,
        x.size,
        endpoint=True,
    )

    def invq(u):
        u = np.asarray(u, dtype=float)
        u = np.clip(
            u,
            1e-12,
            1.0 - 1e-12,
        )
        return np.interp(
            u,
            u_grid,
            x,
        )

    return invq


def qGS_ge_MS_legacy(u, MS):
    """
    Reproduce the earlier conditional-lognormal GS mapping:

        GS ~ fitted lognormal, truncated below at MS

    This is the logic used by the earlier qGS_ge_MS implementation.

    The current canonical model_definitions.py already contains the
    fitted GS lognormal distribution as md._gs_dist.
    """
    u = np.asarray(u, dtype=float)
    MS = np.asarray(MS, dtype=float)

    u = np.clip(
        u,
        1e-12,
        1.0 - 1e-12,
    )

    if not hasattr(md, "_gs_dist"):
        raise AttributeError(
            "model_definitions.py does not expose md._gs_dist. "
            "The diagnostic needs the fitted GS lognormal distribution."
        )

    F0 = md._gs_dist.cdf(MS)
    F0 = np.clip(
        F0,
        0.0,
        1.0 - 1e-12,
    )

    u_adj = (
        F0
        + (1.0 - F0) * u
    )

    u_adj = np.clip(
        u_adj,
        1e-12,
        1.0 - 1e-12,
    )

    return md._gs_dist.ppf(
        u_adj
    )


def evaluate_composite(
    MS,
    GS,
    S,
    CFB,
    CFF,
    legacy_rb_mass=0.0,
):
    """
    Return raw and 100%-bounded composite fnecC.
    """
    MS = np.asarray(MS, dtype=float)
    GS = np.asarray(GS, dtype=float)
    S = np.asarray(S, dtype=float)
    CFB = np.asarray(CFB, dtype=float)
    CFF = np.asarray(CFF, dtype=float)

    if np.any(S <= 0) or np.any(~np.isfinite(S)):
        raise ValueError(
            "All SOC values must be finite and > 0."
        )

    NB = (
        CFB * MS
    )

    if legacy_rb_mass > 0:
        GS_for_NF = np.maximum(
            GS - legacy_rb_mass * MS,
            0.0,
        )
    else:
        GS_for_NF = GS

    NF = (
        CFF * GS_for_NF
    )

    total = (
        NB + NF
    )

    raw = (
        100.0 * total / S
    )

    bounded = np.minimum(
        raw,
        100.0,
    )

    return raw, bounded


def analyze_case(
    problem,
    Y,
    parameter_names,
    seed,
    case,
    gs_mapping,
    response,
):
    """
    Calculate Sobol indices for one diagnostic case.
    """
    if np.any(~np.isfinite(Y)):
        raise RuntimeError(
            f"{case}: non-finite model outputs encountered."
        )

    Si = sobol_analyze(
        problem,
        Y,
        calc_second_order=False,
        num_resamples=100,
        conf_level=0.95,
        print_to_console=False,
        seed=seed,
    )

    return pd.DataFrame({
        "case": case,
        "gs_mapping": gs_mapping,
        "response": response,
        "parameter": parameter_names,
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })


def output_diagnostics(
    case,
    gs_mapping,
    response,
    MS,
    GS,
    S,
    raw,
    Y,
):
    return {
        "case": case,
        "gs_mapping": gs_mapping,
        "response": response,
        "n": len(Y),

        "prop_GS_lt_MS": float(
            np.mean(
                GS < MS
            )
        ),

        "prop_GS_gt_S": float(
            np.mean(
                GS > S
            )
        ),

        "prop_MS_gt_S": float(
            np.mean(
                MS > S
            )
        ),

        "raw_median": float(
            np.median(raw)
        ),

        "raw_p99": float(
            np.quantile(
                raw,
                0.99,
            )
        ),

        "raw_max": float(
            np.max(raw)
        ),

        "prop_raw_gt_100": float(
            np.mean(
                raw > 100.0
            )
        ),

        "target_variance": float(
            np.var(
                Y,
                ddof=1,
            )
        ),

        "target_median": float(
            np.median(Y)
        ),

        "target_max": float(
            np.max(Y)
        ),
    }


# =====================================================================
# Main
# =====================================================================

def main(
    n_base=4096,
    seed=1234,
    lhs_path=None,
    outdir=None,
    legacy_rb_mass=0.0,
):
    if not is_power_of_two(n_base):
        raise ValueError(
            f"n_base={n_base} is not a power of two."
        )

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

    if outdir is None:
        outdir = (
            md.REPO_ROOT
            / "data"
            / "derived"
            / "gsa_diagnostics"
        )
    else:
        outdir = Path(outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not lhs_path.exists():
        raise FileNotFoundError(
            f"Canonical composite LHS file not found: {lhs_path}"
        )

    lhs = pd.read_csv(
        lhs_path
    )

    required = [
        "MS",
        "GS",
        "S",
        "CFB",
        "CFF",
    ]

    missing = [
        x
        for x in required
        if x not in lhs.columns
    ]

    if missing:
        raise ValueError(
            f"{lhs_path} missing required columns: {missing}"
        )

    invq = {
        name: make_invq(
            lhs[name].to_numpy()
        )
        for name in required
    }

    # -------------------------------------------------------------
    # ONE common Sobol design for every case.
    # -------------------------------------------------------------

    parameter_names = [
        "MS",
        "GS",
        "S",
        "CFB",
        "CFF",
    ]

    problem = {
        "num_vars": 5,
        "names": [
            "u_MS",
            "u_GS",
            "u_S",
            "u_CFB",
            "u_CFF",
        ],
        "bounds": [[0.0, 1.0]] * 5,
    }

    U = sobol_sample(
        problem,
        N=n_base,
        calc_second_order=False,
        scramble=True,
        seed=seed,
    )

    expected_rows = (
        n_base * (5 + 2)
    )

    if U.shape != (
        expected_rows,
        5,
    ):
        raise RuntimeError(
            "Unexpected Sobol sample shape. "
            f"Expected {(expected_rows, 5)}, received {U.shape}."
        )

    # -------------------------------------------------------------
    # Shared physical values
    # -------------------------------------------------------------

    MS = invq["MS"](
        U[:, 0]
    )

    S = invq["S"](
        U[:, 2]
    )

    CFB = invq["CFB"](
        U[:, 3]
    )

    CFF = invq["CFF"](
        U[:, 4]
    )

    # Current independent-GS mapping
    GS_independent = invq["GS"](
        U[:, 1]
    )

    # Earlier conditional mapping: guarantees GS >= MS
    GS_ge_MS = qGS_ge_MS_legacy(
        U[:, 1],
        MS,
    )

    # -------------------------------------------------------------
    # Model outputs
    # -------------------------------------------------------------

    raw_ind, bounded_ind = evaluate_composite(
        MS=MS,
        GS=GS_independent,
        S=S,
        CFB=CFB,
        CFF=CFF,
        legacy_rb_mass=legacy_rb_mass,
    )

    raw_ge, bounded_ge = evaluate_composite(
        MS=MS,
        GS=GS_ge_MS,
        S=S,
        CFB=CFB,
        CFF=CFF,
        legacy_rb_mass=legacy_rb_mass,
    )

    cases = [
        (
            "A_independent_GS_raw",
            "independent empirical GS marginal",
            "raw",
            GS_independent,
            raw_ind,
            raw_ind,
        ),
        (
            "B_independent_GS_bounded",
            "independent empirical GS marginal",
            "bounded at 100%",
            GS_independent,
            raw_ind,
            bounded_ind,
        ),
        (
            "C_GS_ge_MS_raw",
            "legacy conditional GS >= MS",
            "raw",
            GS_ge_MS,
            raw_ge,
            raw_ge,
        ),
        (
            "D_GS_ge_MS_bounded",
            "legacy conditional GS >= MS",
            "bounded at 100%",
            GS_ge_MS,
            raw_ge,
            bounded_ge,
        ),
    ]

    sobol_rows = []
    diagnostic_rows = []

    for (
        case,
        gs_mapping,
        response,
        GS_case,
        raw,
        Y,
    ) in cases:

        print()
        print("=" * 72)
        print(case)
        print("=" * 72)

        sobol = analyze_case(
            problem=problem,
            Y=Y,
            parameter_names=parameter_names,
            seed=seed,
            case=case,
            gs_mapping=gs_mapping,
            response=response,
        )

        sobol_rows.append(
            sobol
        )

        diagnostic_rows.append(
            output_diagnostics(
                case=case,
                gs_mapping=gs_mapping,
                response=response,
                MS=MS,
                GS=GS_case,
                S=S,
                raw=raw,
                Y=Y,
            )
        )

        print(
            sobol[
                [
                    "parameter",
                    "S1",
                    "ST",
                ]
            ]
            .sort_values(
                "ST",
                ascending=False,
            )
            .to_string(
                index=False
            )
        )

    sobol_all = pd.concat(
        sobol_rows,
        ignore_index=True,
    )

    diagnostics_all = pd.DataFrame(
        diagnostic_rows
    )

    sobol_path = (
        outdir
        / "composite_four_case_sobol.csv"
    )

    diagnostics_path = (
        outdir
        / "composite_four_case_diagnostics.csv"
    )

    sobol_all.to_csv(
        sobol_path,
        index=False,
    )

    diagnostics_all.to_csv(
        diagnostics_path,
        index=False,
    )

    # -------------------------------------------------------------
    # Compact comparison table for the console
    # -------------------------------------------------------------

    st_compare = (
        sobol_all[
            [
                "case",
                "parameter",
                "ST",
            ]
        ]
        .pivot(
            index="parameter",
            columns="case",
            values="ST",
        )
        .reset_index()
    )

    print()
    print("=" * 72)
    print("TOTAL-ORDER SOBOL INDEX COMPARISON")
    print("=" * 72)
    print(
        st_compare.to_string(
            index=False
        )
    )

    print()
    print("=" * 72)
    print("DIAGNOSTICS")
    print("=" * 72)
    print(
        diagnostics_all[
            [
                "case",
                "prop_GS_lt_MS",
                "prop_GS_gt_S",
                "prop_raw_gt_100",
                "target_variance",
            ]
        ].to_string(
            index=False
        )
    )

    print()
    print("legacy_rb_mass =", legacy_rb_mass)
    print("Saved:")
    print(" ", sobol_path)
    print(" ", diagnostics_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--n_base",
        type=int,
        default=4096,
        help="Sobol base sample size (power of two).",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
    )

    parser.add_argument(
        "--lhs_path",
        default=None,
    )

    parser.add_argument(
        "--outdir",
        default=None,
    )

    parser.add_argument(
        "--legacy_rb_mass",
        type=float,
        default=0.0,
        help=(
            "Optional fixed mass-ratio bacterial GlcN correction. "
            "Use 1.16 to test the older composite formulation."
        ),
    )

    args = parser.parse_args()

    main(
        n_base=args.n_base,
        seed=args.seed,
        lhs_path=args.lhs_path,
        outdir=args.outdir,
        legacy_rb_mass=args.legacy_rb_mass,
    )