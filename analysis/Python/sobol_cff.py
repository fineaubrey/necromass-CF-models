#!/usr/bin/env python3
"""
Sobol GSA for the fungal conversion factor CFF.

Default base sample size: N = 4096.
"""

import argparse
from pathlib import Path

import pandas as pd

from . import model_definitions as md
from . import sobol_common as sc


def transform(U):
    return pd.DataFrame({
        "cF": 300.0 + 400.0 * U[:, 0],
        "GF": md.qGF(U[:, 1]),
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
        d=2,
        n_base=n_base,
        seed=seed,
        calc_second_order=False,
        scramble=True,
    )

    params = transform(U)

    Y = md.CFF_fun(
        params["cF"],
        params["GF"],
    )

    summary, Si = sc.run_sobol_analysis(
        problem=problem,
        y=Y,
        real_names=["cF", "GF"],
        seed=seed,
        num_resamples=num_resamples,
        calc_second_order=False,
    )

    outdir.mkdir(parents=True, exist_ok=True)

    sc.save_sobol_results(
        summary,
        Si,
        outdir,
        "sobol_CFF_saltelli",
    )

    diag = {
        "analysis": "CFF",
        "n_base": n_base,
        "n_model_evaluations": len(Y),
        "seed": seed,
        **sc.summarize_vector(
            Y,
            "CFF",
        ),
    }

    sc.diagnostics_to_frame(
        diag
    ).to_csv(
        outdir / "sobol_CFF_diagnostics.csv",
        index=False,
    )

    if save_draws:
        pd.DataFrame(
            U,
            columns=[
                "u_cF",
                "u_GF",
            ],
        ).to_csv(
            outdir / "sobol_CFF_draws_uniform.csv",
            index=False,
        )

        params.to_csv(
            outdir / "sobol_CFF_draws_params.csv",
            index=False,
        )

        pd.DataFrame({
            "CFF": Y,
        }).to_csv(
            outdir / "sobol_CFF_outputs.csv",
            index=False,
        )

    print("\n[CFF] Sobol indices\n")
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
