#!/usr/bin/env python3
"""
Sobol GSA for fungal necromass carbon NF.

Inputs:
    MS, GS, cF, GF, rB_molar

The bacterial GlcN:MurA ratio is converted from molar to mass units
before it is applied to mass-based soil concentrations.

Default base sample size: N = 4096.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import model_definitions as md
from . import sobol_common as sc


def transform(U):
    rB_molar = (
        4.0 * U[:, 4]
    )

    return pd.DataFrame({
        "MS": md.qMS(U[:, 0]),
        "GS": md.qGS(U[:, 1]),
        "cF": 300.0 + 400.0 * U[:, 2],
        "GF": md.qGF(U[:, 3]),
        "rB_molar": rB_molar,
        "rB_mass": md.rb_molar_to_mass(
            rB_molar
        ),
    })


def run(
    n_base=4096,
    seed=1234,
    outdir=None,
    num_resamples=100,
    save_draws=True,
):
    if outdir is None:
        outdir = md.REPO_ROOT / "data" / "derived" / "gsa"
    else:
        outdir = Path(outdir)

    problem, U = sc.generate_sobol_design(
        d=5,
        n_base=n_base,
        seed=seed,
        calc_second_order=False,
        scramble=True,
    )

    params = transform(U)

    CFF = md.CFF_fun(
        params["cF"],
        params["GF"],
    )

    Y = md.NF_from_mass_ratio(
        CFF=CFF,
        GS=params["GS"],
        MS=params["MS"],
        rB_mass=params["rB_mass"],
    )

    summary, Si = sc.run_sobol_analysis(
        problem=problem,
        y=Y,
        real_names=[
            "MS",
            "GS",
            "cF",
            "GF",
            "rB",
        ],
        seed=seed,
        num_resamples=num_resamples,
        calc_second_order=False,
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    sc.save_sobol_results(
        summary,
        Si,
        outdir,
        "sobol_NF_saltelli",
    )

    diag = {
        "analysis": "NF",
        "n_base": n_base,
        "n_model_evaluations": len(Y),
        "seed": seed,
        "prop_NF_zero": float(
            np.mean(
                np.isclose(
                    Y,
                    0.0,
                )
            )
        ),
        **sc.summarize_vector(
            Y,
            "NF",
        ),
    }

    sc.diagnostics_to_frame(
        diag
    ).to_csv(
        outdir / "sobol_NF_diagnostics.csv",
        index=False,
    )

    if save_draws:
        pd.DataFrame(
            U,
            columns=[
                "u_MS",
                "u_GS",
                "u_cF",
                "u_GF",
                "u_rB",
            ],
        ).to_csv(
            outdir / "sobol_NF_draws_uniform.csv",
            index=False,
        )

        params.to_csv(
            outdir / "sobol_NF_draws_params.csv",
            index=False,
        )

        pd.DataFrame({
            "CFF": CFF,
            "NF": Y,
        }).to_csv(
            outdir / "sobol_NF_outputs.csv",
            index=False,
        )

    print("\n[NF] Sobol indices\n")
    print(
        summary.to_string(
            index=False
        )
    )

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--n_base",
        type=int,
        default=4096,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
    )

    parser.add_argument(
        "--outdir",
        default=None,
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

    run(
        n_base=args.n_base,
        seed=args.seed,
        outdir=args.outdir,
        num_resamples=args.num_resamples,
        save_draws=(
            not args.no_save_draws
        ),
    )
