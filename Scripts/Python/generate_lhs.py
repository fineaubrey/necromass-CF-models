#!/usr/bin/env python3
"""
Generate the canonical LHS uncertainty-propagation datasets.

Repository location
-------------------
scripts/Python/generate_lhs.py

Run from repository root
------------------------
python -m analysis.Python.generate_lhs

or explicitly:

python -m analysis.Python.generate_lhs \
    --n 50000 \
    --seed 1234 \
    --outdir data/derived/lhs

Outputs
-------
data/derived/lhs/
    CFB_lhs.csv
    NF_lhs.csv
    fnecC5_lhs.csv
    fnecC10_lhs.csv

Important
---------
- LHS sample size defaults to N = 50,000.
- Soil MS, GS, and S are sampled from independent fitted marginals.
- Sampled soil states used for fnecC are screened only for the
  molecular-carbon input criterion:
      0.430*MS + 0.402*GS <= S
- No upper bound is imposed on CF-derived fnecC in the LHS analysis.
- rB is used only in the fungal / trait-resolved formulations. It is
  sampled as a MOLAR ratio and converted to a MASS ratio before
  application to mass-based soil concentrations.
- The five-input composite fnecC model applies NO rB correction:
      NB = CFB * MS
      NF = CFF * GS
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import model_definitions as md


# =====================================================================
# 1. LHS helper
# =====================================================================

def lhs_uniform(
    n,
    d,
    seed,
):
    """Simple Latin hypercube sample in the open unit interval."""
    rng = np.random.default_rng(seed)

    cut = np.linspace(
        0.0,
        1.0,
        n + 1,
    )

    u = rng.uniform(
        size=(n, d)
    )

    out = np.empty(
        (n, d),
        dtype=float,
    )

    for j in range(d):

        strata = (
            cut[:-1]
            + u[:, j]
            * (cut[1:] - cut[:-1])
        )

        rng.shuffle(
            strata
        )

        out[:, j] = strata

    return np.clip(
        out,
        1e-12,
        1.0 - 1e-12,
    )


# =====================================================================
# 2. Shared mappings
# =====================================================================

def map_bacterial_traits(U):
    cB = 300.0 + 400.0 * U[:, 0]
    MGP = md.qMGP(U[:, 1])
    MGN = md.qMGN(U[:, 2])
    fGP = 0.1 + 0.8 * U[:, 3]

    return (
        cB,
        MGP,
        MGN,
        fGP,
    )


def map_fungal_traits(U):
    cF = 300.0 + 400.0 * U[:, 0]
    GF = md.qGF(U[:, 1])

    return (
        cF,
        GF,
    )


def draw_feasible_soil_lhs(
    n,
    seed,
    max_rounds=100,
):
    """
    Draw independent MS, GS, and S marginals and retain input states
    satisfying:
        0.430*MS + 0.402*GS <= S

    Additional LHS batches are generated until n feasible soil states
    have been retained.
    """
    accepted = []
    n_accepted = 0

    batch_n = max(
        n,
        int(
            np.ceil(
                1.02 * n
            )
        ),
    )

    for i in range(max_rounds):

        U = lhs_uniform(
            batch_n,
            3,
            seed + i,
        )

        MS = md.qMS(
            U[:, 0]
        )

        GS = md.qGS(
            U[:, 1]
        )

        S = md.qS(
            U[:, 2]
        )

        feasible = md.molecular_carbon_feasible(
            MS,
            GS,
            S,
        )

        accepted.append(
            pd.DataFrame({
                "MS": MS[feasible],
                "GS": GS[feasible],
                "S": S[feasible],
            })
        )

        n_accepted += int(
            feasible.sum()
        )

        if n_accepted >= n:
            break

        remaining = (
            n
            - n_accepted
        )

        batch_n = max(
            remaining,
            int(
                np.ceil(
                    1.10 * remaining
                )
            ),
        )

    if n_accepted < n:
        raise RuntimeError(
            f"Unable to obtain {n} feasible soil LHS draws "
            f"after {max_rounds} batches; obtained {n_accepted}."
        )

    out = pd.concat(
        accepted,
        ignore_index=True,
    ).iloc[:n].copy()

    return out


# =====================================================================
# 3. CFB model
# =====================================================================

def make_cfb_lhs(
    n,
    seed,
):
    U = lhs_uniform(
        n,
        4,
        seed,
    )

    (
        cB,
        MGP,
        MGN,
        fGP,
    ) = map_bacterial_traits(U)

    CFB = md.CFB_fun(
        cB,
        MGP,
        MGN,
        fGP,
    )

    return pd.DataFrame({
        "cB": cB,
        "MGP": MGP,
        "MGN": MGN,
        "fGP": fGP,
        "CFB": CFB,
    })


# =====================================================================
# 4. Fungal necromass model
# =====================================================================

def make_nf_lhs(
    n,
    seed,
):
    U = lhs_uniform(
        n,
        5,
        seed,
    )

    MS = md.qMS(
        U[:, 0]
    )

    GS = md.qGS(
        U[:, 1]
    )

    cF = (
        300.0
        + 400.0 * U[:, 2]
    )

    GF = md.qGF(
        U[:, 3]
    )

    rB_molar = (
        4.0
        * U[:, 4]
    )

    rB_mass = md.rb_molar_to_mass(
        rB_molar
    )

    CFF = md.CFF_fun(
        cF,
        GF,
    )

    NF0 = (
        GS
        * CFF
    )

    NF_rB = md.NF_from_mass_ratio(
        CFF=CFF,
        GS=GS,
        MS=MS,
        rB_mass=rB_mass,
    )

    dNF = (
        NF_rB
        - NF0
    )

    return pd.DataFrame({
        "MS": MS,
        "GS": GS,
        "cF": cF,
        "GF": GF,
        "CFF": CFF,
        "rB_molar": rB_molar,
        "rB_mass": rB_mass,
        "NF0": NF0,
        "NF_rB": NF_rB,
        "dNF": dNF,
    })


# =====================================================================
# 5. Composite five-input fnecC model
#
# IMPORTANT:
#   NO bacterial GlcN correction is applied here.
#   NB = CFB * MS
#   NF = CFF * GS
# =====================================================================

def make_fnecc5_lhs(
    n,
    seed,
):
    soil = draw_feasible_soil_lhs(
        n,
        seed,
    )

    U = lhs_uniform(
        n,
        2,
        seed + 100,
    )

    CFB = md.qCFB(
        U[:, 0]
    )

    CFF = md.qCFF(
        U[:, 1]
    )

    fnecC_raw = md.fnecC_composite_raw(
        MS=soil["MS"].to_numpy(),
        GS=soil["GS"].to_numpy(),
        S=soil["S"].to_numpy(),
        CFB=CFB,
        CFF=CFF,
    )

    return pd.DataFrame({
        "MS": soil["MS"].to_numpy(),
        "GS": soil["GS"].to_numpy(),
        "S": soil["S"].to_numpy(),
        "CFB": CFB,
        "CFF": CFF,
        "input_molecular_C": md.molecular_carbon(
            soil["MS"].to_numpy(),
            soil["GS"].to_numpy(),
        ),
        "fnecC5_pct_uncapped": fnecC_raw,
    })


# =====================================================================
# 6. Trait-resolved ten-input fnecC model
# =====================================================================

def make_fnecc10_lhs(
    n,
    seed,
):
    soil = draw_feasible_soil_lhs(
        n,
        seed,
    )

    U = lhs_uniform(
        n,
        7,
        seed + 200,
    )

    cB = (
        300.0
        + 400.0 * U[:, 0]
    )

    MGP = md.qMGP(
        U[:, 1]
    )

    MGN = md.qMGN(
        U[:, 2]
    )

    fGP = (
        0.1
        + 0.8 * U[:, 3]
    )

    cF = (
        300.0
        + 400.0 * U[:, 4]
    )

    GF = md.qGF(
        U[:, 5]
    )

    rB_molar = (
        4.0
        * U[:, 6]
    )

    rB_mass = md.rb_molar_to_mass(
        rB_molar
    )

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
        soil["MS"].to_numpy(),
    )

    NF = md.NF_from_mass_ratio(
        CFF=CFF,
        GS=soil["GS"].to_numpy(),
        MS=soil["MS"].to_numpy(),
        rB_mass=rB_mass,
    )

    S = soil["S"].to_numpy()

    # Use the canonical shared trait-resolved model definition so this
    # LHS implementation cannot silently diverge from model_definitions.py.
    fnecC_raw = md.fnecC_trait_raw(
        MS=soil["MS"].to_numpy(),
        GS=soil["GS"].to_numpy(),
        S=S,
        cB=cB,
        MGP=MGP,
        MGN=MGN,
        fGP=fGP,
        cF=cF,
        GF=GF,
        rB_molar=rB_molar,
    )

    return pd.DataFrame({
        "MS": soil["MS"].to_numpy(),
        "GS": soil["GS"].to_numpy(),
        "S": S,
        "cB": cB,
        "MGP": MGP,
        "MGN": MGN,
        "fGP": fGP,
        "cF": cF,
        "GF": GF,
        "CFB": CFB,
        "CFF": CFF,
        "rB_molar": rB_molar,
        "rB_mass": rB_mass,
        "NB": NB,
        "NF": NF,
        "input_molecular_C": md.molecular_carbon(
            soil["MS"].to_numpy(),
            soil["GS"].to_numpy(),
        ),
        "fnecC10_pct_uncapped": fnecC_raw,
    })


# =====================================================================
# 7. Main
# =====================================================================

def main(
    n=50_000,
    seed=1234,
    outdir=None,
):
    if outdir is None:
        outdir = (
            md.REPO_ROOT
            / "data"
            / "derived"
            / "lhs"
        )
    else:
        outdir = Path(
            outdir
        )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        f"[LHS] N={n:,}; seed={seed}"
    )

    print(
        "[LHS] output directory:",
        outdir,
    )

    jobs = [
        (
            "CFB_lhs.csv",
            make_cfb_lhs(
                n,
                seed + 1,
            ),
        ),
        (
            "NF_lhs.csv",
            make_nf_lhs(
                n,
                seed + 2,
            ),
        ),
        (
            "fnecC5_lhs.csv",
            make_fnecc5_lhs(
                n,
                seed + 3,
            ),
        ),
        (
            "fnecC10_lhs.csv",
            make_fnecc10_lhs(
                n,
                seed + 4,
            ),
        ),
    ]

    for filename, df in jobs:

        path = (
            outdir
            / filename
        )

        df.to_csv(
            path,
            index=False,
        )

        print(
            "  saved:",
            path,
            f"({len(df):,} rows)",
        )

    print("[LHS] Done.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--n",
        type=int,
        default=50_000,
        help="LHS samples per model (default: 50000).",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
    )

    args = parser.parse_args()

    main(
        n=args.n,
        seed=args.seed,
        outdir=args.outdir,
    )