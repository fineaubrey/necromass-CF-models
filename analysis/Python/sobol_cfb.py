#!/usr/bin/env python3
"""
Sobol GSA for the bacterial conversion factor CFB.

Default base sample size: N = 4096.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import model_definitions as md
from . import sobol_common as sc


def transform(U):
    return pd.DataFrame({
        "cB": 300.0 + 400.0 * U[:, 0],
        "MGP": md.qMGP(U[:, 1]),
        "MGN": md.qMGN(U[:, 2]),
        "fGP": 0.1 + 0.8 * U[:, 3],
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
        d=4,
        n_base=n_base,
        seed=seed,
        calc_second_order=False,
        scramble=True,
    )

    params = transform(U)

    Y = md.CFB_fun(
        params["cB"],
        params["MGP"],
        params["MGN"],
        params["fGP"],
    )

    summary, Si = sc.run_sobol_analysis(
        problem=problem,
        y=Y,
        real_names=["cB", "MGP", "MGN", "fGP"],
        seed=seed,
        num_resamples=num_resamples,
        calc_second_order=False,
    )

    stem = "sobol_CFB_saltelli"
    outdir.mkdir(parents=True, exist_ok=True)

    sc.save_sobol_results(
        summary,
        Si,
        outdir,
        stem,
    )

    diag = {
        "analysis": "CFB",
        "n_base": n_base,
        "n_model_evaluations": len(Y),
        "seed": seed,
        **sc.summarize_vector(
            Y,
            "CFB",
        ),
    }

    sc.diagnostics_to_frame(
        diag
    ).to_csv(
        outdir / "sobol_CFB_diagnostics.csv",
        index=False,
    )

    if save_draws:
        pd.DataFrame(
            U,
            columns=[
                "u_cB",
                "u_MGP",
                "u_MGN",
                "u_fGP",
            ],
        ).to_csv(
            outdir / "sobol_CFB_draws_uniform.csv",
            index=False,
        )

        params.to_csv(
            outdir / "sobol_CFB_draws_params.csv",
            index=False,
        )

        pd.DataFrame({
            "CFB": Y,
        }).to_csv(
            outdir / "sobol_CFB_outputs.csv",
            index=False,
        )

    print("\n[CFB] Sobol indices\n")
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
