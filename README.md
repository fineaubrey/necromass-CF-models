# Necromass Conversion Factor Models

Reproducible analysis code for the manuscript:

**“Conversion factors as a source of uncertainty in amino sugar-based estimates of soil microbial necromass carbon”**

This repository evaluates how microbial trait variability, conversion-factor parameterization, and soil amino-sugar measurements influence estimates of microbial necromass carbon in soils.

The workflow has three main components:

1. **Microbial trait variability** — quantify variation in bacterial muramic acid (MurA) and fungal glucosamine (GlcN) traits.
2. **Uncertainty propagation and global sensitivity analysis** — propagate microbial-trait uncertainty and soil-input variability through conversion-factor models using Latin hypercube sampling (LHS) and Sobol sensitivity analysis.
3. **Global-soil application** — compare conversion-factor estimates with molecule-level accounting of carbon directly contained in measured amino sugars across a global soil dataset.

---

## Repository structure

```text
necromass-CF-models/
├── README.md
├── .gitignore
│
├── data/
│   ├── processed/
│   │   ├── Dataset1.csv
│   │   └── Dataset2.csv
│   │
│   └── derived/
│       ├── lhs/
│       ├── gsa/
│       └── global_soils/
│
├── analysis/
│   ├── diagnostics/
│   └── figures/
│
├── Scripts/
│   ├── Python/
│   │   ├── __init__.py
│   │   ├── archive/
│   │   │   ├── sobol_fnecC5_saltelli.py
│   │   │   └── sobol_fnecC10_saltelli.py
│   │   ├── model_definitions.py
│   │   ├── generate_lhs.py
│   │   ├── summarize_lhs.py
│   │   ├── sobol_common.py
│   │   ├── sobol_cfb.py
│   │   ├── sobol_cff.py
│   │   ├── sobol_nf.py
│   │   ├── sobol_fnecC5_saltelli.py
│   │   ├── sobol_fnecC10_saltelli.py
│   │   └── run_gsa.py
│   │
│   ├── R/
│   │   ├── 01_traits/
│   │   │   ├── 00_setup.R
│   │   │   ├── 01_trait_summary.R
│   │   │   ├── 02_trait_data.R
│   │   │   ├── 02_trait_statistics.R
│   │   │   └── 03_Fig2_trait_distributions.R
│   │   │
│   │   ├── 02_uncertainty_sensitivity/
│   │   │   └── 03_Fig3_sobol_indices.R
│   │   │
│   │   └── 03_global_soils/
│   │       ├── 01_prepare_global_soils.R
│   │       ├── 02_Fig4_global_soils_model_comparisons.R
│   │       ├── 03_model_comparison_statistics.R
│   │       ├── 04_agreement_CCC_bland_altman.R
│   │       └── 05_Fig5_structural_feasibility.R
│
├── results/
│   └── tables/
│
└── manuscript/
    └── figures/
```

Generated intermediate files are written to `data/derived/`. Manuscript-ready figures are written to `manuscript/figures/`, and tabular results are written to `results/tables/`.

---

## Software requirements

### R

The analysis uses R packages including:

```r
tidyverse
here
ggplot2
scales
patchwork
boot
DescTools
viridis
```

Additional packages may be required by individual trait-analysis scripts.

### Python

The Python workflow requires:

```text
numpy
pandas
scipy
SALib
```

Run all Python commands from the repository root.

---

## Input datasets

The workflow expects two processed input files:

```text
data/processed/Dataset1.csv
data/processed/Dataset2.csv
```

### Dataset 1

Dataset 1 contains microbial trait observations used to parameterize bacterial and fungal conversion factors.

The canonical trait dataset includes:

- 444 bacterial MurA observations
- 781 fungal GlcN observations
- Oomycota excluded from the fungal dataset

### Dataset 2

Dataset 2 contains soil MurA, GlcN, SOC, and associated metadata used in the empirical global-soil analysis.

In the source file:

- MurA is stored in `ug g^-1 soil`
- GlcN is stored in `ug g^-1 soil`
- SOC is stored in `mg C g^-1 soil`

`01_prepare_global_soils.R` converts MurA and GlcN to `mg g^-1 soil` once and saves the canonical analytical dataset. Downstream scripts should **not divide MurA or GlcN by 1000 again**.

---

# Analysis workflow

## 1. Microbial trait variability

Open the repository as an RStudio project and run:

```r
source("Scripts/R/01_traits/01_trait_summary.R")
source("Scripts/R/01_traits/02_trait_statistics.R")
source("Scripts/R/01_traits/03_Fig2_trait_distributions.R")
```

These scripts summarize trait distributions, test taxonomic differences, and generate Figure 2.

---

## 2. Latin hypercube uncertainty propagation

From PowerShell, the RStudio Terminal, or another shell opened at the repository root:

```bash
python -m Scripts.Python.generate_lhs
```

Then summarize the LHS outputs:

```bash
python -m Scripts.Python.summarize_lhs
```

Canonical LHS outputs are written to:

```text
data/derived/lhs/
```

including:

```text
CFB_lhs.csv
NF_lhs.csv
fnecC5_lhs.csv
fnecC10_lhs.csv
LHS_summary_statistics.csv
```

### Final LHS choices

The canonical uncertainty-propagation workflow uses:

- `N = 50,000`
- independent fitted soil marginal distributions
- MurA and SOC truncated to observed minimum-maximum ranges
- GlcN sampled over its selected empirical range
- a molecular-carbon feasibility condition:

```math
0.430M_S + 0.402G_S \leq S
```

No `GS >= MS` constraint is imposed.

The LHS uncertainty-propagation analysis does **not** cap CF-derived `fnecC` at 100% SOC.

---

## 3. Sobol global sensitivity analysis

Run all canonical Sobol analyses with:

```bash
python -m Scripts.Python.run_gsa
```

Default Sobol base sample sizes are:

| Analysis | Base N |
|---|---:|
| Bacterial conversion factor (`CFB`) | 4096 |
| Fungal conversion factor (`CFF`) | 4096 |
| Fungal necromass (`NF`) | 4096 |
| Composite `fnecC` | 4096 |
| Trait-resolved `fnecC` | 8192 |

Outputs are written to:

```text
data/derived/gsa/
```

### Final Sobol model choices

For the **composite five-input model**:

```math
N_B = CF_B M_S
```

```math
N_F = CF_F G_S
```

```math
f_{necC} = 100\frac{N_B + N_F}{S}
```

- `MS` and `GS` are varied independently.
- No `GS >= MS` constraint is imposed.
- No bacterial GlcN (`rB`) correction is used.

For the **trait-resolved ten-input model**:

```math
N_B = CF_BM_S
```

```math
G_{S,\mathrm{corr}} =
\max(G_S-r_{B,\mathrm{mass}}M_S,0)
```

```math
N_F = CF_FG_{S,\mathrm{corr}}
```

where the bacterial GlcN:MurA ratio is sampled on a molar basis and converted before it is applied to mass-based soil concentrations:

```math
r_{B,\mathrm{mass}}
=
r_{B,\mathrm{molar}}
\frac{MW_{\mathrm{GlcN}}}{MW_{\mathrm{MurA}}}
```

For the SOC-normalized Sobol analyses, all Sobol rows are retained and the primary sensitivity target is the raw, unbounded response:

```math
f_{necC,\mathrm{raw}}
=
100\frac{N_B+N_F}{S}
```

As a robustness analysis, the sensitivity analysis is repeated after imposing a physical upper bound:

```math
f_{necC,\mathrm{bounded}}
=
\min(f_{necC,\mathrm{raw}},100)
```

Raw, unbounded results are used in the main-text analysis and Figure 3. Results after bounding at 100% of SOC are reported in Appendix D and used to evaluate whether imposing a physical upper limit alters the principal parameter rankings.

### Generate Figure 3

After running the GSA:

```r
source("Scripts/R/02_uncertainty_sensitivity/03_Fig3_sobol_indices.R")
```

---

## 4. Global-soil empirical analysis

### Step 1 — Prepare the analytical dataset

```r
source("Scripts/R/03_global_soils/01_prepare_global_soils.R")
```

This script:

- converts MurA and GlcN from `ug g^-1` to `mg g^-1`
- retains observations with finite MurA, GlcN, and SOC
- requires `SOC > 0`
- does not impose `GS >= MS`
- does not apply statistical outlier filtering
- saves the canonical analytical sample to:

```text
data/derived/global_soils/Dataset2_analytical_sample.csv
```

### Step 2 — Apply empirical models and generate Figure 4

```r
source("Scripts/R/03_global_soils/02_Fig4_global_soils_model_comparisons.R")
```

Four empirical quantities are calculated:

#### Standard conversion-factor model

```math
CF_{2006}
=
100\frac{45M_S+9G_S}{S}
```

#### Hu et al. (2024) formulation

```math
CF_{2024}
=
100
\frac{
31.3M_S+
10.8\max(G_S-1.16M_S,0)
}{S}
```

The fixed `1.16` bacterial GlcN correction is retained here because it is part of the published CF2024 formulation being compared. It is **not** used in the five-input composite GSA.

#### Direct MurA + GlcN molecular carbon

```math
C_{\mathrm{MurA+GlcN}}
=
100
\frac{
0.430M_S+0.402G_S
}{S}
```

#### Direct GlcN molecular carbon

```math
C_{\mathrm{GlcN}}
=
100
\frac{
0.402G_S
}{S}
```

Empirical model outputs are **not capped at 100% SOC**. Values above 100% are retained and reported as diagnostics.

Canonical model estimates are saved to:

```text
data/derived/global_soils/global_soils_model_estimates_wide.csv
data/derived/global_soils/global_soils_model_estimates_long.csv
```

---

## 5. Paired numerical comparison among empirical outputs

Run:

```r
source("Scripts/R/03_global_soils/03_model_comparison_statistics.R")
```

The script summarizes prespecified paired numerical differences among model outputs. It reports median paired differences with 10,000-resample bootstrap 95% confidence intervals. These comparisons quantify the magnitude and direction of numerical differences; they are not treated as tests of whether the models estimate the same underlying quantity.

Results are written to:

```text
results/tables/
```

The exact output filenames are defined in `03_model_comparison_statistics.R`.

---

## 6. Agreement analysis

Run:

```r
source("Scripts/R/03_global_soils/04_agreement_CCC_bland_altman.R")
```

Agreement diagnostics are restricted to comparisons between formulations intended to represent comparable quantities. CF-based necromass estimates and molecule-level carbon calculations are not interpreted as interchangeable estimates of the same pool.

The script calculates:

- Lin's concordance correlation coefficient (CCC)
- 10,000-resample bootstrap CCC confidence intervals
- Bland-Altman mean difference
- 95% limits of agreement
- bootstrap confidence intervals for the mean difference and limits of agreement

For Bland-Altman comparisons:

```math
\mathrm{difference}
=
\mathrm{method1}-\mathrm{method2}
```

so a positive mean difference means method 1 tends to produce larger estimates.

Outputs include:

```text
results/tables/global_soils_CCC_pairwise.csv
results/tables/global_soils_bland_altman_summary.csv
data/derived/global_soils/global_soils_bland_altman_long.csv
analysis/figures/BlandAltman_grid.png
manuscript/figures/BlandAltman_grid.pdf
```

---

## 7. Structural-feasibility figure

Run:

```r
source("Scripts/R/03_global_soils/05_Fig5_structural_feasibility.R")
```

This figure evaluates the relationship among soil GlcN, SOC, and CF-derived necromass estimates.

Panel A shows reference iso-`fnecC` contours for the fungal scaling term:

```math
f_{necC}
=
100\frac{CF_FG_S}{S}
```

using `CFF = 9`.

These contours are a structural reference and are **not** the exact full CF2006 surface because the full CF2006 model also contains the bacterial term `45 * MurA`.

Panel B shows the full empirical CF2006 estimate as a function of soil GlcN.

---

# Reproducing the analysis

A complete run can be performed in the following order.

### R: traits

```r
source("Scripts/R/01_traits/01_trait_summary.R")
source("Scripts/R/01_traits/02_trait_statistics.R")
source("Scripts/R/01_traits/03_Fig2_trait_distributions.R")
```

### Python: LHS and GSA

```bash
python -m Scripts.Python.generate_lhs
python -m Scripts.Python.summarize_lhs
python -m Scripts.Python.run_gsa
```

### R: sensitivity figure

```r
source("Scripts/R/02_uncertainty_sensitivity/03_Fig3_sobol_indices.R")
```

### R: global-soil analyses

```r
source("Scripts/R/03_global_soils/01_prepare_global_soils.R")
source("Scripts/R/03_global_soils/02_Fig4_global_soils_model_comparisons.R")
source("Scripts/R/03_global_soils/03_model_comparison_statistics.R")
source("Scripts/R/03_global_soils/04_agreement_CCC_bland_altman.R")
source("Scripts/R/03_global_soils/05_Fig5_structural_feasibility.R")
```

---

## Reproducibility notes

Several model definitions were explored during development. The scripts in the main analysis workflow implement the final choices used for the reproducible analysis.

In particular:

- Oomycota are excluded from the fungal trait dataset.
- `GS >= MS` is **not** imposed as a sampling constraint.
- The five-input composite model contains **no `rB` correction**.
- The trait-resolved model does include an `rB` correction after molar-to-mass conversion.
- Classical Sobol coordinates are varied independently through empirical marginal distributions.
- The primary SOC-normalized Sobol analysis uses raw, unbounded outputs; bounding at 100% is evaluated as a robustness analysis.
- Empirical global-soil estimates are not capped at 100%.
- Dataset 2 amino-sugar concentrations are converted from `ug g^-1` to `mg g^-1` only once, during analytical-sample preparation.
- Diagnostic and legacy scripts should be kept separate from the canonical workflow.

---

## Citation

Please cite the associated manuscript when using this code or derived outputs.
