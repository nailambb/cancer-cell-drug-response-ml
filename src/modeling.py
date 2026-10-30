import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_regression
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, GroupKFold, KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

SEED = 42
ROOT = Path(__file__).resolve().parent.parent
REPORTS = ROOT / "reports"


def get_cv(scheme, n_splits=5):
    """'kfold' = random split (optimistic); 'group' = whole lineages held out (harder)."""
    if scheme == "kfold":
        return KFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    if scheme == "group":
        return GroupKFold(n_splits=n_splits)
    raise ValueError(scheme)


def build_pipeline(model_name):
    """Everything that learns from data lives inside the pipeline, so it is
    refit on each training fold and never sees the test fold."""
    models = {
        "dummy": DummyRegressor(strategy="mean"),
        "ridge": Ridge(),
        "elasticnet": ElasticNet(max_iter=5000, random_state=SEED),
        "rf": RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1),
    }
    return Pipeline([
        ("var", VarianceThreshold(threshold=0.01)),
        ("select", SelectKBest(f_regression, k=1000)),
        ("scale", StandardScaler()),
        ("model", models[model_name]),
    ])


def get_param_grid(model_name):
    """Keys are '<step>__<param>'. Keep grids small: nested CV multiplies the cost."""
    k = {"select__k": [500, 2000]}
    grids = {
        "dummy": {},
        "ridge": {**k, "model__alpha": [10, 100, 1000, 10000]},
        "elasticnet": {**k, "model__alpha": [0.01, 0.1, 1],
                       "model__l1_ratio": [0.2, 0.5, 0.8]},
        "rf": {"select__k": [500], "model__max_features": ["sqrt", 0.3]},
    }
    return grids[model_name]

def evaluate_model(model_name, X, y, groups, scheme):
    """Nested CV: the outer split evaluates, GridSearchCV inside tunes on the
    training fold only. Returns metrics plus out-of-fold predictions."""
    inner = KFold(n_splits=3, shuffle=True, random_state=SEED)
    tuned = GridSearchCV(build_pipeline(model_name), get_param_grid(model_name),
                         cv=inner, scoring="neg_root_mean_squared_error", n_jobs=-1)
    # groups only goes to the OUTER splitter
    preds = cross_val_predict(tuned, X, y, cv=get_cv(scheme), groups=groups)

    metrics = {
        "model": model_name,
        "scheme": scheme,
        "r2": r2_score(y, preds),
        "rmse": float(np.sqrt(mean_squared_error(y, preds))),
        "spearman": spearmanr(y, preds).correlation,
    }
    return metrics, pd.Series(preds, index=y.index, name=f"{model_name}_{scheme}")

def lineage_baseline(y, groups, scheme):
    """Predict each held-out line as the mean AUC of its lineage in the training folds."""
    preds = pd.Series(index=y.index, dtype=float, name=f"lineage_mean_{scheme}")
    for tr, te in get_cv(scheme).split(y, y, groups):
        means = y.iloc[tr].groupby(groups.iloc[tr]).mean()
        preds.iloc[te] = groups.iloc[te].map(means).fillna(y.iloc[tr].mean()).values
    metrics = {
        "model": "lineage_mean",
        "scheme": scheme,
        "r2": r2_score(y, preds),
        "rmse": float(np.sqrt(mean_squared_error(y, preds))),
        "spearman": spearmanr(y, preds).correlation,
    }
    return metrics, preds


def run_all(drug_id, drug_name, X, y, groups,
            model_names=("dummy", "ridge", "elasticnet", "rf")):
    REPORTS.mkdir(exist_ok=True)
    rows, oof = [], {}
    for scheme in ("kfold", "group"):
        for name in model_names:
            print(f"{drug_name} | {scheme} | {name}")
            metrics, preds = evaluate_model(name, X, y, groups, scheme)
            rows.append(metrics)
            oof[preds.name] = preds
        # lineage-only baseline, same splits, same output format
        metrics, preds = lineage_baseline(y, groups, scheme)
        rows.append(metrics)
        oof[preds.name] = preds
        
    results = pd.DataFrame(rows)
    results.insert(0, "drug", drug_name)
    results.to_csv(REPORTS / f"results_{drug_name}.csv", index=False)
    pd.DataFrame(oof).assign(y=y, lineage=groups).to_csv(
        REPORTS / f"oof_{drug_name}.csv")
    return results


def get_feature_importance(model_name, X, y):
    """Refit tuned pipeline on ALL data (for interpretation only, not evaluation).
    Returns a Series of coefficients/importances indexed by gene name."""
    tuned = GridSearchCV(build_pipeline(model_name), get_param_grid(model_name),
                         cv=KFold(3, shuffle=True, random_state=SEED),
                         scoring="neg_root_mean_squared_error", n_jobs=-1).fit(X, y)
    best = tuned.best_estimator_
    # Track which genes survive each selection step, or gene names get mislabeled
    names = X.columns[best.named_steps["var"].get_support()]
    names = names[best.named_steps["select"].get_support()]
    model = best.named_steps["model"]
    values = model.coef_ if hasattr(model, "coef_") else model.feature_importances_
    return pd.Series(values, index=names).sort_values(key=abs, ascending=False)


if __name__ == "__main__":
    from data_preprocessing import build_dataset

    drugs = {1373: "Dabrafenib", 1080: "Paclitaxel"}
    drug_id = int(sys.argv[1]) if len(sys.argv) > 1 else 1373
    X, y, groups = build_dataset(drug_id)
    print(run_all(drug_id, drugs[drug_id], X, y, groups))