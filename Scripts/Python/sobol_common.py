"""
Shared utilities for Sobol global sensitivity analyses.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from SALib.analyze.sobol import analyze as sobol_analyze
from SALib.sample.sobol import sample as sobol_sample


def is_power_of_two(n):
    return n > 0 and (n & (n - 1) == 0)


def make_invq(empirical_values, min_n=50):
    """
    Build a linear empirical inverse-CDF function.
    """
    x = np.asarray(empirical_values, dtype=float)
    x = x[np.isfinite(x)]
    x = np.sort(x)

    if x.size < min_n:
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


def generate_sobol_design(
    d,
    n_base,
    seed,
    calc_second_order=False,
    scramble=True,
):
    """
    Generate a SALib Sobol/Saltelli design in an independent unit hypercube.
    """
    if not is_power_of_two(n_base):
        raise ValueError(
            f"n_base={n_base} is not a power of two."
        )

    problem = {
        "num_vars": d,
        "names": [
            f"u{i + 1}"
            for i in range(d)
        ],
        "bounds": [
            [0.0, 1.0]
        ] * d,
    }

    U = sobol_sample(
        problem,
        N=n_base,
        calc_second_order=calc_second_order,
        scramble=scramble,
        seed=seed,
    )

    expected_rows = (
        n_base * (2 * d + 2)
        if calc_second_order
        else n_base * (d + 2)
    )

    if U.shape != (
        expected_rows,
        d,
    ):
        raise RuntimeError(
            "Unexpected Sobol sample shape. "
            f"Expected {(expected_rows, d)}, received {U.shape}."
        )

    return problem, U


def run_sobol_analysis(
    problem,
    y,
    real_names,
    seed,
    num_resamples=100,
    calc_second_order=False,
):
    """
    Calculate first- and total-order Sobol indices.
    """
    y = np.asarray(
        y,
        dtype=float,
    )

    if np.any(
        ~np.isfinite(y)
    ):
        raise ValueError(
            "Sobol target contains non-finite values."
        )

    Si = sobol_analyze(
        problem,
        y,
        calc_second_order=calc_second_order,
        num_resamples=num_resamples,
        conf_level=0.95,
        print_to_console=False,
        seed=seed,
    )

    summary = pd.DataFrame({
        "param": real_names,
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })

    return summary, Si


def save_sobol_results(
    summary,
    Si,
    outdir,
    stem,
):
    """
    Save Sobol summary CSV and raw SALib result arrays.
    """
    outdir = Path(outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary_path = (
        outdir
        / f"{stem}_summary.csv"
    )

    npz_path = (
        outdir
        / f"{stem}_Si.npz"
    )

    summary.to_csv(
        summary_path,
        index=False,
    )

    np.savez(
        npz_path,
        **Si,
    )

    return (
        summary_path,
        npz_path,
    )


def summarize_vector(
    x,
    prefix,
):
    x = np.asarray(
        x,
        dtype=float,
    )

    x = x[
        np.isfinite(x)
    ]

    return {
        f"{prefix}_min": np.min(x),
        f"{prefix}_mean": np.mean(x),
        f"{prefix}_median": np.median(x),
        f"{prefix}_sd": np.std(
            x,
            ddof=1,
        ),
        f"{prefix}_p90": np.quantile(
            x,
            0.90,
        ),
        f"{prefix}_p95": np.quantile(
            x,
            0.95,
        ),
        f"{prefix}_p99": np.quantile(
            x,
            0.99,
        ),
        f"{prefix}_p995": np.quantile(
            x,
            0.995,
        ),
        f"{prefix}_max": np.max(x),
    }


def diagnostics_to_frame(diagnostics):
    return pd.DataFrame({
        "metric": list(
            diagnostics.keys()
        ),
        "value": list(
            diagnostics.values()
        ),
    })
