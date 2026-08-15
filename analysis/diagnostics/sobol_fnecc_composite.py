#!/usr/bin/env python3
"""
Sobol GSA for the five-input composite fnecC model.

Inputs:
    MS, GS, S, CFB, CFF

Marginal distributions are taken from the canonical LHS output:
    data/derived/lhs/fnecC5_lhs.csv

The Sobol design varies those five marginals independently.

Raw model:
    fnecC_raw = 100 * (CFB*MS + CFF*GS) / S

Sobol target:
    fnecC_bounded = min(fnecC_raw, 100)

No Sobol rows are removed.

Default base sample size: N = 4096.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import model_definitions as md
from . import sobol_common as sc


def run(
    n_base=4096,
    seed=1234,
    outdir=None,
    lhs_path=None,
    num_resamples=100,
    save_draws=True,
):
    if outdir is None:
        outdir = md.REPO_ROOT / "data" / "derived" / "gsa"
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

    lhs = pd.read_csv(
        lhs_path
    )

    names = [
        "MS",
        "GS",
        "S",
        "CFB",
        "CFF",
    ]

    missing = [
        x
        for x in names
        if x not in lhs.columns
    ]

    if missing:
        raise ValueError(
            f"{lhs_path} missing columns: {missing}"
        )

    invq = {
        name: sc.make_invq(
            lhs[name].to_numpy()
        )
        for name in names
    }

    problem, U = sc.generate_sobol_design(
        d=5,
        n_base=n_base,
        seed=seed,
        calc_second_order=False,
        scramble=True,
    )

    params = pd.DataFrame({
        name: invq[name](
            U[:, i]
        )
        for i, name in enumerate(
            names
        )
    })

    raw = md.fnecC_composite_raw(
        MS=params["MS"],
        GS=params["GS"],
        S=params["S"],
        CFB=params["CFB"],
        CFF=params["CFF"],
    )

    bounded = np.minimum(
        raw,
        100.0,
    )

    summary, Si = sc.run_sobol_analysis(
        problem=problem,
        y=bounded,
        real_names=names,
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
        "sobol_fnecC5_saltelli",
    )

    diag = {
        "analysis": "fnecC composite",
        "n_base": n_base,
        "n_model_evaluations": len(raw),
        "seed": seed,
        "marginal_source": str(
            lhs_path
        ),
        "prop_raw_gt_100": float(
            np.mean(
                raw > 100.0
            )
        ),
        "prop_at_upper_bound": float(
            np.mean(
                np.isclose(
                    bounded,
                    100.0,
                )
            )
        ),
        **sc.summarize_vector(
            raw,
            "raw",
        ),
        **sc.summarize_vector(
            bounded,
            "bounded",
        ),
    }

    sc.diagnostics_to_frame(
        diag
    ).to_csv(
        outdir / "sobol_fnecC5_diagnostics.csv",
        index=False,
    )

    if save_draws:
        pd.DataFrame(
            U,
            columns=[
                "u_MS",
                "u_GS",
                "u_S",
                "u_CFB",
                "u_CFF",
            ],
        ).to_csv(
            outdir / "sobol_fnecC5_draws_uniform.csv",
            index=False,
        )

        params.to_csv(
            outdir / "sobol_fnecC5_draws_params.csv",
            index=False,
        )

        pd.DataFrame({
            "fnecC5_raw": raw,
            "fnecC5_bounded": bounded,
        }).to_csv(
            outdir / "sobol_fnecC5_outputs.csv",
            index=False,
        )

    print("\n[fnecC composite] Sobol indices\n")
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
        "--lhs_path",
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
        lhs_path=args.lhs_path,
        num_resamples=args.num_resamples,
        save_draws=(
            not args.no_save_draws
        ),
    )
