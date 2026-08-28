"""
Canonical model definitions for the necromass-CF uncertainty analysis.

Repository location
-------------------
scripts/Python/model_definitions.py

Canonical inputs
----------------
data/processed/Dataset1.csv
data/processed/Dataset2.csv

Key analysis choices
--------------------
- Dataset 1 fungal trait distribution excludes Oomycota (N = 781).
- MGP, MGN, and GF are sampled from empirical P10-P90 quantile ranges.
- MS and GS marginals are restricted to their empirical P10-P90 ranges.
- S is truncated to its observed range.
- Soil marginals are varied independently in the current LHS workflow.
- rB is used ONLY in the trait-resolved formulation. It is sampled as
  a MOLAR GlcN:MurA ratio and converted to a MASS ratio before application
  to mass-based soil concentrations.
- Molecular-carbon feasibility uses exact molecular carbon fractions:
      MurA = 0.430
      GlcN = 0.402
- LHS fnecC outputs are left uncapped.
- Sobol scripts may apply the physical output bound min(NB + NF, S).

The composite fnecC workflow applies no explicit bacterial GlcN
correction, matching Eq. 3a and its accompanying Methods text.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


# =====================================================================
# 0. Repository paths and constants
# =====================================================================

REPO_ROOT = Path(__file__).resolve().parents[2]

P_TRAITS = REPO_ROOT / "data" / "processed" / "Dataset1.csv"
P_SOILS = REPO_ROOT / "data" / "processed" / "Dataset2.csv"

EPS = 1e-12

TRAIT_QLO = 0.10
TRAIT_QHI = 0.90

# Soil biomarker sampling bounds
SOIL_MS_QLO = 0.10
SOIL_MS_QHI = 0.90

SOIL_GS_QLO = 0.10
SOIL_GS_QHI = 0.90

# Molecular weights (g mol^-1)
MW_GLCN = 179.17
MW_MURA = 251.23
RB_MOLAR_TO_MASS = MW_GLCN / MW_MURA

# Exact carbon mass fractions used in the current analysis.
MURA_C_FRACTION = 0.430
GLCN_C_FRACTION = 0.402


# =====================================================================
# 1. Robust CSV loading
# =====================================================================

def read_csv_robust(path):
    """Read a CSV using a small sequence of plausible encodings."""
    path = Path(path)

    for enc in ("utf-8", "utf-8-sig", "cp1252", "latin1"):
        try:
            return pd.read_csv(path, encoding=enc)
        except UnicodeDecodeError:
            continue

    return pd.read_csv(
        path,
        encoding="latin1",
        encoding_errors="replace",
    )


def _require_columns(df, required, label):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            f"{label} is missing required columns: {missing}"
        )


# =====================================================================
# 2. Dataset 1 — canonical microbial trait vectors
# =====================================================================

if not P_TRAITS.exists():
    raise FileNotFoundError(
        f"Dataset1.csv not found at: {P_TRAITS}"
    )

traits_raw = read_csv_robust(P_TRAITS)
traits_raw = traits_raw.replace([np.inf, -np.inf], np.nan)

_require_columns(
    traits_raw,
    ["Kingdom", "Phylum", "Gram", "MurA", "GlcN"],
    "Dataset1",
)

# Remove the stray encoding artifact previously found in a phylum label.
traits_raw["Phylum"] = (
    traits_raw["Phylum"]
    .astype("string")
    .str.replace("Â", "", regex=False)
)

kingdom = (
    traits_raw["Kingdom"]
    .astype("string")
    .str.strip()
    .str.lower()
)

phylum = (
    traits_raw["Phylum"]
    .astype("string")
    .str.strip()
    .str.lower()
)

gram = (
    traits_raw["Gram"]
    .astype("string")
    .str.strip()
    .str.upper()
)

is_bacteria = kingdom.eq("bacteria")

# Canonical fungal definition from Phase 1:
# exact Fungi rows, excluding Oomycota.
is_fungi = kingdom.eq("fungi") & ~phylum.eq("oomycota")

murA_bacteria = pd.to_numeric(
    traits_raw.loc[is_bacteria, "MurA"],
    errors="coerce",
).dropna().to_numpy(dtype=float)

murA_GP_vec = pd.to_numeric(
    traits_raw.loc[is_bacteria & gram.eq("GP"), "MurA"],
    errors="coerce",
).dropna().to_numpy(dtype=float)

murA_GN_vec = pd.to_numeric(
    traits_raw.loc[is_bacteria & gram.eq("GN"), "MurA"],
    errors="coerce",
).dropna().to_numpy(dtype=float)

glcN_F_vec = pd.to_numeric(
    traits_raw.loc[is_fungi, "GlcN"],
    errors="coerce",
).dropna().to_numpy(dtype=float)

if murA_bacteria.size != 444:
    raise ValueError(
        "Canonical bacterial MurA dataset should contain 444 "
        f"observations; found {murA_bacteria.size}."
    )

if glcN_F_vec.size != 781:
    raise ValueError(
        "Canonical fungal GlcN dataset should contain 781 "
        "observations after excluding Oomycota; "
        f"found {glcN_F_vec.size}."
    )

if murA_GP_vec.size != 262 or murA_GN_vec.size != 182:
    raise ValueError(
        "Unexpected Gram-group counts. "
        f"GP={murA_GP_vec.size}, GN={murA_GN_vec.size}."
    )


# =====================================================================
# 3. Dataset 2 — canonical soil variables
# =====================================================================

if not P_SOILS.exists():
    raise FileNotFoundError(
        f"Dataset2.csv not found at: {P_SOILS}"
    )

soil_raw = read_csv_robust(P_SOILS)
soil_raw = soil_raw.replace([np.inf, -np.inf], np.nan)

_require_columns(
    soil_raw,
    ["MurA", "GlcN", "SOC"],
    "Dataset2",
)

soil_numeric = pd.DataFrame({
    # Dataset2 stores amino sugars in micrograms g^-1 soil.
    "MS": pd.to_numeric(
        soil_raw["MurA"],
        errors="coerce",
    ) / 1000.0,

    "GS": pd.to_numeric(
        soil_raw["GlcN"],
        errors="coerce",
    ) / 1000.0,

    # SOC is already mg C g^-1 soil.
    "S": pd.to_numeric(
        soil_raw["SOC"],
        errors="coerce",
    ),
})

soils = (
    soil_numeric
    .replace([np.inf, -np.inf], np.nan)
    .dropna()
)

# Fitted distributions require strictly positive values.
soils = soils[
    (soils["MS"] > 0)
    & (soils["GS"] > 0)
    & (soils["S"] > 0)
].copy()

if soils.empty:
    raise ValueError(
        "No positive, finite soil observations remain after cleaning."
    )


# =====================================================================
# 4. Empirical trait quantile functions
# =====================================================================

def q_emp(u, arr):
    """Empirical inverse-CDF mapping."""
    u = np.asarray(u, dtype=float)
    u = np.clip(u, 1e-12, 1.0 - 1e-12)

    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]

    return np.quantile(arr, u)


def q_emp_bounded(
    u,
    arr,
    qlo=TRAIT_QLO,
    qhi=TRAIT_QHI,
):
    """
    Map u in [0,1] through the empirical distribution restricted
    to the requested quantile interval without discarding observations.
    """
    u = np.asarray(u, dtype=float)
    u = np.clip(u, 0.0, 1.0)

    u_bounded = qlo + (qhi - qlo) * u

    return q_emp(
        u_bounded,
        arr,
    )


def qMGP(u):
    return q_emp_bounded(
        u,
        murA_GP_vec,
    )


def qMGN(u):
    return q_emp_bounded(
        u,
        murA_GN_vec,
    )


def qGF(u):
    return q_emp_bounded(
        u,
        glcN_F_vec,
    )


# =====================================================================
# 5. Fitted soil marginal distributions
# =====================================================================

def fit_lognorm_params(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x) & (x > 0)]

    log_x = np.log(x)

    return (
        float(log_x.mean()),
        float(log_x.std(ddof=1)),
    )


ms_meanlog, ms_sdlog = fit_lognorm_params(
    soils["MS"].to_numpy()
)

gs_meanlog, gs_sdlog = fit_lognorm_params(
    soils["GS"].to_numpy()
)

S_shape, _, S_scale = stats.gamma.fit(
    soils["S"].to_numpy(dtype=float),
    floc=0,
)

_ms_dist = stats.lognorm(
    s=ms_sdlog,
    scale=np.exp(ms_meanlog),
)

_gs_dist = stats.lognorm(
    s=gs_sdlog,
    scale=np.exp(gs_meanlog),
)

_s_dist = stats.gamma(
    a=S_shape,
    loc=0,
    scale=S_scale,
)

# Current LHS domain:
# MS and GS -> empirical P10-P90
# S -> observed min-max

MS_LO, MS_HI = np.quantile(
    soils["MS"].to_numpy(dtype=float),
    [SOIL_MS_QLO, SOIL_MS_QHI],
)

MS_LO = float(MS_LO)
MS_HI = float(MS_HI)

GS_LO, GS_HI = np.quantile(
    soils["GS"].to_numpy(dtype=float),
    [SOIL_GS_QLO, SOIL_GS_QHI],
)

GS_LO = float(GS_LO)
GS_HI = float(GS_HI)

S_LO = float(soils["S"].min())
S_HI = float(soils["S"].max())


def q_truncated(u, dist, lo, hi):
    """Inverse-CDF mapping for a fitted distribution truncated to [lo, hi]."""
    u = np.asarray(u, dtype=float)
    u = np.clip(u, 1e-12, 1.0 - 1e-12)

    f_lo = dist.cdf(lo)
    f_hi = dist.cdf(hi)

    u_adj = f_lo + u * (f_hi - f_lo)

    return dist.ppf(u_adj)


def qMS(u):
    return q_truncated(
        u,
        _ms_dist,
        MS_LO,
        MS_HI,
    )


def qGS(u):
    return q_truncated(
        u,
        _gs_dist,
        GS_LO,
        GS_HI,
    )


def qS(u):
    return q_truncated(
        u,
        _s_dist,
        S_LO,
        S_HI,
    )


# =====================================================================
# 6. Ratio conversion and physical-input diagnostics
# =====================================================================

def rb_molar_to_mass(rB_molar):
    """Convert GlcN:MurA molar ratio to the mass ratio applied to soils."""
    return (
        np.asarray(rB_molar, dtype=float)
        * RB_MOLAR_TO_MASS
    )


def molecular_carbon(
    MS,
    GS,
):
    """Carbon directly contained in measured MurA + GlcN."""
    MS = np.asarray(MS, dtype=float)
    GS = np.asarray(GS, dtype=float)

    return (
        MURA_C_FRACTION * MS
        + GLCN_C_FRACTION * GS
    )


def molecular_carbon_feasible(
    MS,
    GS,
    S,
):
    """
    Input-state feasibility criterion used in LHS:
    0.430*MS + 0.402*GS <= S.
    """
    S = np.asarray(S, dtype=float)

    return molecular_carbon(
        MS,
        GS,
    ) <= S


# =====================================================================
# 7. Conversion-factor and necromass models
# =====================================================================

def CFB_fun(
    cB,
    MGP,
    MGN,
    fGP,
):
    denominator = (
        np.asarray(MGP, dtype=float)
        * np.asarray(fGP, dtype=float)
        + np.asarray(MGN, dtype=float)
        * (1.0 - np.asarray(fGP, dtype=float))
    )

    denominator = np.maximum(
        denominator,
        EPS,
    )

    return (
        np.asarray(cB, dtype=float)
        / denominator
    )


def CFF_fun(
    cF,
    GF,
):
    GF = np.maximum(
        np.asarray(GF, dtype=float),
        EPS,
    )

    return (
        np.asarray(cF, dtype=float)
        / GF
    )


def NB_fun(
    CFB,
    MS,
):
    return (
        np.asarray(CFB, dtype=float)
        * np.asarray(MS, dtype=float)
    )


def NF_from_mass_ratio(
    CFF,
    GS,
    MS,
    rB_mass,
):
    corrected_GS = np.maximum(
        np.asarray(GS, dtype=float)
        - np.asarray(rB_mass, dtype=float)
        * np.asarray(MS, dtype=float),
        0.0,
    )

    return (
        corrected_GS
        * np.asarray(CFF, dtype=float)
    )


def NF_from_molar_ratio(
    CFF,
    GS,
    MS,
    rB_molar,
):
    return NF_from_mass_ratio(
        CFF=CFF,
        GS=GS,
        MS=MS,
        rB_mass=rb_molar_to_mass(rB_molar),
    )


# =====================================================================
# 8. SOC-normalized model outputs
# =====================================================================

def fnecC_composite_raw(
    MS,
    GS,
    S,
    CFB,
    CFF,
):
    """
    Five-input composite fnecC implementation.

    The uncertain inputs are MS, GS, S, CFB, and CFF.

    No explicit bacterial GlcN correction is applied in the composite
    formulation. Thus:

        NB = CFB * MS
        NF = CFF * GS
        fnecC = 100 * (NB + NF) / S

    This matches Eq. 3a and the accompanying Methods text.
    """
    S_safe = np.maximum(
        np.asarray(S, dtype=float),
        EPS,
    )

    NB = NB_fun(
        CFB,
        MS,
    )

    NF = (
        np.asarray(CFF, dtype=float)
        * np.asarray(GS, dtype=float)
    )

    total_necromass = NB + NF

    return (
        100.0
        * total_necromass
        / S_safe
    )


def fnecC_trait_raw(
    MS,
    GS,
    S,
    cB,
    MGP,
    MGN,
    fGP,
    cF,
    GF,
    rB_molar,
):
    """Ten-input trait-resolved fnecC implementation."""
    CFB = CFB_fun(
        cB,
        MGP,
        MGN,
        fGP,
    )

    CFF = CFF_fun(
        cF,
        GF,
    )

    NB = NB_fun(
        CFB,
        MS,
    )

    NF = NF_from_molar_ratio(
        CFF=CFF,
        GS=GS,
        MS=MS,
        rB_molar=rB_molar,
    )

    S_safe = np.maximum(
        np.asarray(S, dtype=float),
        EPS,
    )

    return (
        100.0
        * (NB + NF)
        / S_safe
    )


def bound_total_necromass(
    NB,
    NF,
    S,
):
    """
    Physical output bound used only for SOC-normalized Sobol analyses:
        N_total,bounded = min(NB + NF, S)
    """
    S_safe = np.maximum(
        np.asarray(S, dtype=float),
        EPS,
    )

    return np.minimum(
        np.asarray(NB, dtype=float)
        + np.asarray(NF, dtype=float),
        S_safe,
    )


# =====================================================================
# 9. Empirical composite-CF distributions
# =====================================================================

def build_empirical_CFs(
    n=200_000,
    seed=123,
):
    """
    Propagate the underlying trait distributions to obtain empirical
    composite CFB and CFF pools. qCFB/qCFF then sample their P10-P90
    ranges for the five-input composite model.
    """
    rng = np.random.default_rng(seed)

    u = rng.random(
        (n, 6)
    )

    cB = 300.0 + 400.0 * u[:, 0]
    MGP = qMGP(u[:, 1])
    MGN = qMGN(u[:, 2])
    fGP = 0.1 + 0.8 * u[:, 3]

    cF = 300.0 + 400.0 * u[:, 4]
    GF = qGF(u[:, 5])

    CFB = CFB_fun(
        cB,
        MGP,
        MGN,
        fGP,
    )

    CFF = CFF_fun(
        cF,
        GF,
    )

    CFB = CFB[
        np.isfinite(CFB)
        & (CFB > 0)
    ]

    CFF = CFF[
        np.isfinite(CFF)
        & (CFF > 0)
    ]

    return (
        np.sort(CFB),
        np.sort(CFF),
    )


CFB_emp, CFF_emp = build_empirical_CFs()


def _q_cf_bounded(
    u,
    pool,
    qlo=0.10,
    qhi=0.90,
):
    u = np.asarray(u, dtype=float)
    u = np.clip(u, 0.0, 1.0)

    u_bounded = (
        qlo
        + (qhi - qlo) * u
    )

    return np.quantile(
        pool,
        u_bounded,
    )


def qCFB(u):
    return _q_cf_bounded(
        u,
        CFB_emp,
    )


def qCFF(u):
    return _q_cf_bounded(
        u,
        CFF_emp,
    )


# =====================================================================
# 10. Lightweight diagnostics
# =====================================================================

if __name__ == "__main__":

    print("Repository root:", REPO_ROOT)
    print()

    print("Canonical trait counts:")
    print("  bacterial MurA:", murA_bacteria.size)
    print("  GP MurA:", murA_GP_vec.size)
    print("  GN MurA:", murA_GN_vec.size)
    print("  fungal GlcN:", glcN_F_vec.size)
    print()

    print("Trait P10-P90:")
    print(
        "  MGP:",
        np.quantile(
            murA_GP_vec,
            [0.10, 0.90],
        ),
    )
    print(
        "  MGN:",
        np.quantile(
            murA_GN_vec,
            [0.10, 0.90],
        ),
    )
    print(
        "  GF:",
        np.quantile(
            glcN_F_vec,
            [0.10, 0.90],
        ),
    )
    print()

    print("Soil sampling domain:")
    print("  MS P10-P90:", MS_LO, MS_HI)
    print("  GS P10-P90:", GS_LO, GS_HI)
    print("  S observed range:", S_LO, S_HI)