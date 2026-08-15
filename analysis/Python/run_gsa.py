#!/usr/bin/env python3
"""
Run the complete canonical Sobol global sensitivity workflow.

Run from repository root
------------------------
python -m analysis.Python.run_gsa

Canonical analyses
------------------
1. Bacterial conversion factor (CFB)
2. Fungal conversion factor (CFF)
3. Fungal necromass (NF)
4. Five-input composite fnecC
5. Ten-input trait-resolved fnecC

Final methodological choices
----------------------------
- CFB, CFF, and NF use the canonical standalone Sobol scripts.
- Composite fnecC:
    * MS and GS sampled independently
    * NO GS >= MS constraint
    * NO rB correction
    * bounded fnecC is the Sobol target
- Trait-resolved fnecC:
    * MS and GS sampled independently
    * NO GS >= MS constraint
    * rB correction IS applied
    * rB sampled as a molar ratio and converted to a mass ratio
    * bounded fnecC is the Sobol target
- No structured Sobol rows are deleted.

Default base sample sizes
-------------------------
CFB             4096
CFF             4096
NF              4096
Composite       4096
Trait-resolved  8192

Outputs
-------
data/derived/gsa/

Expected summary files:
    sobol_CFB_saltelli_summary.csv
    sobol_CFF_saltelli_summary.csv
    sobol_NF_saltelli_summary.csv
    sobol_fnecC5_saltelli_summary.csv
    sobol_fnecC10_saltelli_summary.csv
"""

import argparse
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------
# Repository paths
# ---------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]

GSA_DIR = (
    REPO_ROOT
    / "data"
    / "derived"
    / "gsa"
)


# ---------------------------------------------------------------------
# Canonical modules and expected outputs
# ---------------------------------------------------------------------

ANALYSES = [
    {
        "label": "CFB",
        "module": "analysis.Python.sobol_cfb",
        "default_n": 4096,
        "summary": "sobol_CFB_saltelli_summary.csv",
    },
    {
        "label": "CFF",
        "module": "analysis.Python.sobol_cff",
        "default_n": 4096,
        "summary": "sobol_CFF_saltelli_summary.csv",
    },
    {
        "label": "NF",
        "module": "analysis.Python.sobol_nf",
        "default_n": 4096,
        "summary": "sobol_NF_saltelli_summary.csv",
    },
    {
        "label": "Composite fnecC",
        "module": "analysis.Python.sobol_fnecC5_saltelli",
        "default_n": 4096,
        "summary": "sobol_fnecC5_saltelli_summary.csv",
    },
    {
        "label": "Trait-resolved fnecC",
        "module": "analysis.Python.sobol_fnecC10_saltelli",
        "default_n": 8192,
        "summary": "sobol_fnecC10_saltelli_summary.csv",
    },
]


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def run_module(module, n_base):
    """
    Run one canonical Sobol module using the current Python interpreter.
    """
    cmd = [
        sys.executable,
        "-m",
        module,
        "--n_base",
        str(n_base),
    ]

    print()
    print("=" * 78)
    print("RUNNING:", " ".join(cmd))
    print("=" * 78)
    print()

    subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        check=True,
    )


def check_expected_outputs():
    """
    Verify that each canonical summary file exists after the run.
    """
    print()
    print("=" * 78)
    print("OUTPUT CHECK")
    print("=" * 78)

    missing = []

    for analysis in ANALYSES:
        path = (
            GSA_DIR
            / analysis["summary"]
        )

        if path.exists():
            print(
                "[OK] ",
                analysis["label"],
                ": ",
                path.relative_to(REPO_ROOT),
                sep="",
            )
        else:
            print(
                "[MISSING] ",
                analysis["label"],
                ": ",
                path.relative_to(REPO_ROOT),
                sep="",
            )
            missing.append(path)

    if missing:
        raise FileNotFoundError(
            "One or more expected GSA summary files were not created."
        )


# ---------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------

def main(
    n_base_standard=4096,
    n_base_trait=8192,
    skip_cfb=False,
    skip_cff=False,
    skip_nf=False,
    skip_composite=False,
    skip_trait=False,
):
    """
    Run the canonical Sobol GSA workflow.
    """

    GSA_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    skip_lookup = {
        "CFB": skip_cfb,
        "CFF": skip_cff,
        "NF": skip_nf,
        "Composite fnecC": skip_composite,
        "Trait-resolved fnecC": skip_trait,
    }

    print("=" * 78)
    print("CANONICAL GLOBAL SENSITIVITY ANALYSIS")
    print("=" * 78)
    print("Repository root :", REPO_ROOT)
    print("Output directory:", GSA_DIR)
    print("Standard base N :", n_base_standard)
    print("Trait base N    :", n_base_trait)
    print()
    print("Final fnecC choices:")
    print("  Composite:      independent GS/MS; NO rB correction")
    print("  Trait-resolved: independent GS/MS; rB correction YES")
    print("  Sobol target:   fnecC bounded at 100%")
    print()

    for analysis in ANALYSES:

        label = analysis["label"]

        if skip_lookup[label]:
            print(
                "[SKIP]",
                label,
            )
            continue

        n_base = (
            n_base_trait
            if label == "Trait-resolved fnecC"
            else n_base_standard
        )

        run_module(
            module=analysis["module"],
            n_base=n_base,
        )

    # Only require outputs for analyses that were not explicitly skipped.
    print()
    print("=" * 78)
    print("OUTPUT CHECK")
    print("=" * 78)

    missing = []

    for analysis in ANALYSES:

        label = analysis["label"]

        if skip_lookup[label]:
            continue

        path = (
            GSA_DIR
            / analysis["summary"]
        )

        if path.exists():
            print(
                "[OK] ",
                label,
                ": ",
                path.relative_to(REPO_ROOT),
                sep="",
            )
        else:
            print(
                "[MISSING] ",
                label,
                ": ",
                path.relative_to(REPO_ROOT),
                sep="",
            )
            missing.append(path)

    if missing:
        raise FileNotFoundError(
            "One or more expected GSA summary files were not created."
        )

    print()
    print("=" * 78)
    print("GSA WORKFLOW COMPLETE")
    print("=" * 78)
    print()
    print(
        "Next step: regenerate Figure 3 from the summary CSV files in"
    )
    print(
        "  data/derived/gsa/"
    )


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Run all canonical Sobol global sensitivity analyses."
        )
    )

    parser.add_argument(
        "--n_base_standard",
        type=int,
        default=4096,
        help=(
            "Base sample size for CFB, CFF, NF, and composite fnecC "
            "(default: 4096)."
        ),
    )

    parser.add_argument(
        "--n_base_trait",
        type=int,
        default=8192,
        help=(
            "Base sample size for trait-resolved fnecC "
            "(default: 8192)."
        ),
    )

    parser.add_argument(
        "--skip_cfb",
        action="store_true",
    )

    parser.add_argument(
        "--skip_cff",
        action="store_true",
    )

    parser.add_argument(
        "--skip_nf",
        action="store_true",
    )

    parser.add_argument(
        "--skip_composite",
        action="store_true",
    )

    parser.add_argument(
        "--skip_trait",
        action="store_true",
    )

    args = parser.parse_args()

    main(
        n_base_standard=args.n_base_standard,
        n_base_trait=args.n_base_trait,
        skip_cfb=args.skip_cfb,
        skip_cff=args.skip_cff,
        skip_nf=args.skip_nf,
        skip_composite=args.skip_composite,
        skip_trait=args.skip_trait,
    )