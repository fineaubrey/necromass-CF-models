#!/usr/bin/env python3
"""
Summarize canonical LHS uncertainty-propagation outputs.

Run from repository root
------------------------
python -m analysis.Python.summarize_lhs

Output
------
data/derived/lhs/LHS_summary_statistics.csv
"""

from pathlib import Path

import pandas as pd

from . import model_definitions as md


def summarize(x):
    x = pd.Series(x).dropna()

    return {
        "N": len(x),
        "Min": x.min(),
        "Q1": x.quantile(0.25),
        "Median": x.median(),
        "Q3": x.quantile(0.75),
        "Max": x.max(),
        "Mean": x.mean(),
        "SD": x.std(ddof=1),
    }


def main(
    indir=None,
):
    if indir is None:
        indir = (
            md.REPO_ROOT
            / "data"
            / "derived"
            / "lhs"
        )
    else:
        indir = Path(
            indir
        )

    targets = {
        "CFB": {
            "file": "CFB_lhs.csv",
            "column": "CFB",
            "units": "g C g MurA^-1",
        },
        "CFF": {
            "file": "NF_lhs.csv",
            "column": "CFF",
            "units": "g C g GlcN^-1",
        },
        "NF": {
            "file": "NF_lhs.csv",
            "column": "NF_rB",
            "units": "mg C g soil^-1",
        },
        "fnecC composite": {
            "file": "fnecC5_lhs.csv",
            "column": "fnecC5_pct_uncapped",
            "units": "% SOC",
        },
        "fnecC trait-resolved": {
            "file": "fnecC10_lhs.csv",
            "column": "fnecC10_pct_uncapped",
            "units": "% SOC",
        },
    }

    rows = []

    for model, info in targets.items():

        path = (
            indir
            / info["file"]
        )

        if not path.exists():
            raise FileNotFoundError(
                f"Missing LHS output: {path}"
            )

        df = pd.read_csv(
            path
        )

        if info["column"] not in df.columns:
            raise ValueError(
                f"{info['column']} not found in {path}."
            )

        stats = summarize(
            df[
                info["column"]
            ]
        )

        stats["Model"] = model
        stats["Units"] = info["units"]

        rows.append(
            stats
        )

    summary = pd.DataFrame(
        rows
    )[
        [
            "Model",
            "Units",
            "N",
            "Min",
            "Q1",
            "Median",
            "Q3",
            "Max",
            "Mean",
            "SD",
        ]
    ]

    outpath = (
        indir
        / "LHS_summary_statistics.csv"
    )

    summary.to_csv(
        outpath,
        index=False,
    )

    print(
        summary.round(3).to_string(
            index=False
        )
    )

    print(
        "\nSaved:",
        outpath,
    )


if __name__ == "__main__":
    main()
