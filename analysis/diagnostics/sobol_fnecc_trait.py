#!/usr/bin/env python3
"""
Sobol GSA for the ten-input trait-resolved fnecC model.

Inputs:
    MS, GS, S, cB, MGP, MGN, fGP, cF, GF, rB_molar

Marginal distributions are taken from:
    data/derived/lhs/fnecC10_lhs.csv

The Sobol design varies the ten marginals independently.

The rB input remains a molar ratio and is converted to a mass ratio
inside the model evaluation.

Raw model:
    fnecC_raw = 100 * (NB + NF) / S

Sobol target:
    fnecC_bounded = min(fnecC_raw, 100)

No Sobol rows are removed.

Default base sample size: N = 8192.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import model_definitions as md
from . import sobol_common as sc


def run(
    n_base=8192,
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
            / "fnecC10_lhs.csv"
        )
    else:
        lhs_path = Path(lhs_path)

    lhs = pd.read_csv(
        lhs_path
    )

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
        x
        for x in source_names
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
        for name in source_names
    }

    problem, U = sc.generate_sobol_design(
        d=10,
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
            source_names
        )
    })

    CFB = md.CFB_fun(
        params["cB"],
        params["MGP"],
        params["MGN"],
        params["fGP"],
    )

    CFF = md.CFF_fun(
        params["cF"],
        params["GF"],
    )

    NB = md.NB_fun(
        CFB,
        params["MS"],
    )

    rB_mass = md.rb_molar_to_mass(
        params["rB_molar"]
    )

    NF = md.NF_from_mass_ratio(
        CFF=CFF,
        GS=params["GS"],
        MS=params["MS"],
        rB_mass=rB_mass,
    )

    raw = (
        100.0
        * (NB + NF)
        / np.maximum(
            params["S"].to_numpy(
                dtype=float
            ),
            md.EPS,
        )
    )

    bounded = np.minimum(
        raw,
        100.0,
    )

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

    summary, Si = sc.run_sobol_analysis(
        problem=problem,
        y=bounded,
        real_names=real_names,
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
        "sobol_fnecC10_saltelli",
    )

    diag = {
        "analysis": "fnecC trait-resolved",
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
        outdir / "sobol_fnecC10_diagnostics.csv",
        index=False,
    )

    if save_draws:
        pd.DataFrame(
            U,
            columns=[
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
            ],
        ).to_csv(
            outdir / "sobol_fnecC10_draws_uniform.csv",
            index=False,
        )

        params_out = params.copy()
        params_out["rB_mass"] = (
            rB_mass
        )
        params_out["CFB"] = CFB
        params_out["CFF"] = CFF

        params_out.to_csv(
            outdir / "sobol_fnecC10_draws_params.csv",
            index=False,
        )

        pd.DataFrame({
            "NB": NB,
            "NF": NF,
            "fnecC10_raw": raw,
            "fnecC10_bounded": bounded,
        }).to_csv(
            outdir / "sobol_fnecC10_outputs.csv",
            index=False,
        )

    print("\n[fnecC trait-resolved] Sobol indices\n")
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
        default=8192,
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
