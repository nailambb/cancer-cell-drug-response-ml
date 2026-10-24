# Cancer Cell Line Drug Response Prediction

Predicting cancer cell line drug sensitivity (AUC) from baseline gene expression
using GDSC2 drug response data and DepMap expression data.

## Overview

**Question:** can baseline RNA-seq expression predict how sensitive a cancer cell
line is to a given drug, and does that hold up across tissue types?

**Approach:** regularized linear models and tree-based models, evaluated with
cross-validation, including tissue-grouped CV to test generalization to unseen
tissues. Model coefficients and SHAP values are checked against known drug biology.

## Data

Raw data is **not included** in this repo. Download the files below and place
them in `data/raw/`.

| File | Source | Version |
|---|---|---|
| `GDSC2_fitted_dose_response_27Oct23.xlsx` | [GDSC downloads](https://cog.sanger.ac.uk/cmp/download) (GDSC2 IC50 data) | 27Oct23 |
| `Model.csv` | [DepMap portal](https://depmap.org/portal/download) | DepMap Public 26Q3 |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | [DepMap portal](https://depmap.org/portal/download) | DepMap Public 26Q3 |

- **Labels:** `AUC` from GDSC2 (lower = more sensitive). AUC is used over IC50
  because IC50 is censored at the maximum tested concentration.
- **Features:** log2(TPM+1) expression of protein-coding genes (already
  log-transformed; not transformed again).
- **ID mapping:** GDSC `SANGER_MODEL_ID` → DepMap `SangerModelID` → `ModelID`
  → expression row.

### Data processing notes

- One expression profile per model, using `IsDefaultEntryForModel == "Yes"`.
- **Cell line funnel:** 969 GDSC2 cell lines → 966 mapped to `Model.csv`
  → 722 with expression data.
- `SIDM00117` maps to two DepMap models (NCI-H322, NCI-H322M). GDSC screened
  NCI-H322M, which has no expression, so this cell line was excluded rather than
  risk a label/feature mismatch.
- Drugs are identified by `DRUG_ID`, not name: some names (e.g. Docetaxel) appear
  under multiple IDs with different screening conditions.
- No duplicate (drug, cell line) rows were found.
- TODO: note how name mismatches (~3% of mapped lines) were handled.

## Drug selection

See `notebooks/01_drug_selection_eda.ipynb`. All candidate drugs had ~700 cell
lines and 27 lineages, so selection used AUC variance, response distribution shape,
and known biology.

| Drug (ID) | Type | Role | Notes |
|---|---|---|---|
| Dabrafenib (1373) | Targeted (BRAF) | Primary | Mostly resistant with a sensitive tail; expected biomarker BRAF V600E |
| Paclitaxel (1080) | Cytotoxic | Primary | No single driver; tests limits of expression-based prediction |
| Afatinib (1032) | Targeted (EGFR) | Optional extension | EGFR pathway |

Rejected: Gefitinib (almost all lines resistant, AUC std 0.045), Docetaxel 1007
(narrow response range), Docetaxel 1819 (different assay condition, fewer lines).

## Methods

TODO: preprocessing (variance filtering), models (Ridge, Lasso, Elastic Net,
Random Forest, XGBoost), CV strategy (KFold and tissue-grouped GroupKFold, with
lineages under 10 lines merged into "Other"), metrics (R², RMSE, Spearman ρ).
Feature selection and scaling are fit inside CV folds to avoid leakage.

## Results

TODO: results table and figures from `reports/figures/`.

## Limitations

TODO: e.g. cell lines are not patients; expression-only features; small sensitive
subgroups for targeted drugs; models may partly learn tissue of origin.

## Repo structure

```
data/raw/            # downloaded files (gitignored)
data/processed/      # generated outputs (gitignored)
notebooks/           # 01 drug selection EDA, 02 results exploration
src/                 # data_preprocessing.py, modeling.py
reports/figures/     # saved figures
requirements.txt
```

## Reproducing

```bash
pip install -r requirements.txt
# download the three files above into data/raw/
python src/modeling.py
```