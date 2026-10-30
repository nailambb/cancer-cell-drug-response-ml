from functools import lru_cache
from pathlib import Path

import pandas as pd

# Resolve paths from this file so it works from the repo root or from notebooks/
ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

GDSC_FILE = RAW / "GDSC2_fitted_dose_response_27Oct23.xlsx"
MODEL_FILE = RAW / "Model.csv"
EXPR_FILE = RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"

# SIDM00117 maps to NCI-H322 and NCI-H322M. GDSC screened H322M, which has no
# expression, so no valid label/feature pair exists -> exclude both.
EXCLUDE_MODEL_IDS = ["ACH-000837", "ACH-002172"]
MIN_LINEAGE_SIZE = 10

EXPR_META_COLS = ["Unnamed: 0", "SequencingID", "ModelConditionID",
                  "IsDefaultEntryForMC", "IsDefaultEntryForModel"]


@lru_cache(maxsize=1)
def load_gdsc():
    """GDSC2 fitted dose-response. The Excel read is slow, so cache as parquet."""
    cache = PROCESSED / "gdsc2.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    df = pd.read_excel(GDSC_FILE)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache)
    return df


def load_model_metadata():
    """ModelID <-> SangerModelID bridge plus tissue lineage."""
    meta = pd.read_csv(
        MODEL_FILE,
        usecols=["ModelID", "SangerModelID", "CellLineName", "OncotreeLineage"],
    )
    meta = meta.dropna(subset=["SangerModelID"])
    meta = meta[~meta["ModelID"].isin(EXCLUDE_MODEL_IDS)]
    assert meta["SangerModelID"].is_unique, "Duplicate Sanger IDs remain in Model.csv"
    return meta


@lru_cache(maxsize=1)
def load_expression():
    """Full expression matrix: one row per model (default profile), indexed by ModelID."""
    expr = pd.read_csv(EXPR_FILE)
    expr = expr[expr["IsDefaultEntryForModel"] == "Yes"]
    assert expr["ModelID"].is_unique, "ModelID not unique after default filter"
    expr = expr.set_index("ModelID").drop(columns=EXPR_META_COLS, errors="ignore")
    expr = expr.astype("float32")   # halves memory vs float64
    print("Expression matrix:", expr.shape, "| NaN cells:", int(expr.isna().sum().sum()))
    return expr


def harmonize_ids(response, meta, expr):
    """Attach ModelID and lineage to GDSC rows; keep only models with expression."""
    n0 = response["SANGER_MODEL_ID"].nunique()
    merged = response.merge(
        meta,
        left_on="SANGER_MODEL_ID",
        right_on="SangerModelID",
        how="inner",
        validate="many_to_one",
    )
    n1 = merged["ModelID"].nunique()
    merged = merged[merged["ModelID"].isin(expr.index)]
    n2 = merged["ModelID"].nunique()
    print(f"Funnel: {n0} GDSC lines -> {n1} mapped -> {n2} with expression")
    return merged


def merge_rare_lineages(lineage, min_size=MIN_LINEAGE_SIZE):
    """Lineages with fewer than min_size lines become 'Other' (for grouped CV)."""
    counts = lineage.value_counts()
    rare = counts[counts < min_size].index
    return lineage.where(~lineage.isin(rare), "Other")


def build_dataset(drug_id):
    """Return X (models x genes), y (AUC), groups (lineage) for one DRUG_ID."""
    gdsc = load_gdsc()
    resp = gdsc[gdsc["DRUG_ID"] == drug_id].dropna(subset=["AUC"])
    merged = harmonize_ids(resp, load_model_metadata(), load_expression())

    assert merged["ModelID"].is_unique, "Duplicate (drug, cell line) rows"
    merged = merged.set_index("ModelID")

    X = load_expression().loc[merged.index]
    y = merged["AUC"].astype(float)
    groups = merge_rare_lineages(merged["OncotreeLineage"])

    assert (X.index == y.index).all() and (X.index == groups.index).all()
    return X, y, groups


if __name__ == "__main__":
    X, y, groups = build_dataset(1373)   # Dabrafenib
    print(X.shape, y.shape)
    print(groups.value_counts())
