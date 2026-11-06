"""Fit leakage-conscious RNA models for the selected drug panel.

Feature selection is restricted to training data and performance is compared with matched mean and cancer-category baselines.

Run from the repository root with ``python 03_rna_models.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# RNA prediction and matched baselines
# ----------------------------------------------------------------------------

import json
import os
import sys
import time
import zipfile
import warnings
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from scipy.stats import spearmanr, pearsonr

try:
    from xgboost import XGBRegressor

    HAS_XGB = True
except ImportError:
    HAS_XGB = False

# ---------------- USER CONFIGURATION ----------------
ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase2"
SEED = 20261028
MAX_DRUGS = 15  # pilot benchmark: expand after reviewing real results
MIN_TRAIN_PER_DRUG = 120  # excludes drugs with sparse training observations
N_GENES = 256  # highest variability, selected from TRAINING lines only
XGB_ESTIMATORS = 160
RF_ESTIMATORS = 100
CHUNK_GENES = 900  # streamed RNA matrix, prevents oversized temporary arrays
# ----------------------------------------------------

OUT.mkdir(parents=True, exist_ok=True)
REPORT = []
ISSUES = []
START = time.time()


def log(section, msg):
    REPORT.append(f"[{section}] {msg}")


def issue(section, msg):
    ISSUES.append(f"{section}: {msg}")
    log("WARNING", f"{section}: {msg}")


def get_file(prefix):
    files = sorted(
        p
        for p in DATA.iterdir()
        if p.is_file() and p.name.lower().startswith(prefix.lower())
    )
    if len(files) != 1:
        raise FileNotFoundError(
            f"Expected exactly one file starting {prefix!r}; found {[p.name for p in files]}"
        )
    return files[0]


def norm(s):
    x = s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    return x.mask(x.isin(["", "<NA>", "NAN", "NONE"]))


def metrics(y, pred):
    y, pred = np.asarray(y, float), np.asarray(pred, float)
    mask = np.isfinite(y) & np.isfinite(pred)
    y, pred = y[mask], pred[mask]
    if len(y) < 3:
        return {
            "n": len(y),
            "MAE": np.nan,
            "RMSE": np.nan,
            "R2": np.nan,
            "Pearson": np.nan,
            "Spearman": np.nan,
        }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = float(pearsonr(y, pred).statistic) if np.std(y) and np.std(pred) else np.nan
        s = (
            float(spearmanr(y, pred).statistic)
            if np.std(y) and np.std(pred)
            else np.nan
        )
    return {
        "n": len(y),
        "MAE": float(mean_absolute_error(y, pred)),
        "RMSE": float(np.sqrt(mean_squared_error(y, pred))),
        "R2": float(r2_score(y, pred)),
        "Pearson": p,
        "Spearman": s,
    }


def format_selected_results(selected, cancer_baseline):
    """Return the established Phase 2 result schema without changing metric values."""
    result = selected.copy()
    result["selected_model"] = result["model"]
    result["baseline_MAE"] = result["drug_id"].map(cancer_baseline["MAE"])
    result["MAE_improvement"] = result["baseline_MAE"] - result["MAE"]
    result["improved"] = result["MAE_improvement"] > 0
    return result


def pick_group_splits(models):
    # Retain different derived cell lines from the same patient in ONE partition.
    patient = (
        models["patient_id"].astype("string").str.strip()
        if "patient_id" in models
        else pd.Series(pd.NA, index=models.index, dtype="string")
    )
    patient = patient.mask(patient.isin(["", "nan", "None", "<NA>"]))
    group = pd.Series(
        np.where(patient.notna(), "P:" + patient.fillna(""), "M:" + models["model_id"]),
        index=models.index,
    )
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=SEED)
    trainval_idx, test_idx = next(splitter.split(models, groups=group))
    sub = models.iloc[trainval_idx]
    # 0.15/(0.70+0.15) of remaining samples -> approximately 70/15/15.
    splitter2 = GroupShuffleSplit(
        n_splits=1, test_size=(0.15 / 0.85), random_state=SEED + 1
    )
    tr_inner, va_inner = next(splitter2.split(sub, groups=group.iloc[trainval_idx]))
    train_idx, val_idx = trainval_idx[tr_inner], trainval_idx[va_inner]
    names = np.full(len(models), "test", dtype=object)
    names[train_idx] = "train"
    names[val_idx] = "val"
    d = models[["model_id", "patient_id", "tissue", "cancer_type"]].copy()
    d["group_key"] = group.to_numpy()
    d["split"] = names
    assert d.groupby("group_key")["split"].nunique().eq(1).all()
    return d


def read_expression(tpm_zip, matched_ids, training_ids):
    """Read specific wide-file columns and gene rows; metadata rows are not genes."""
    with zipfile.ZipFile(tpm_zip) as z:
        members = [
            m
            for m in z.namelist()
            if m.lower().endswith(".csv") and "rsem_tpm_" in m.lower()
        ]
        if len(members) != 1:
            raise ValueError(f"Expected 1 RSEM TPM matrix; found {members}")
        member = members[0]
        with z.open(member) as stream:
            header = pd.read_csv(stream, nrows=0).columns.to_list()
    sidm = [c for c in header if c.startswith("SIDM") and c in matched_ids]
    if len(sidm) < 100:
        raise ValueError(f"RNA/GDSC shared model count unexpectedly small: {len(sidm)}")
    # The archive has one row-label column; model metadata occupy the first
    # three rows, not three identifier columns.
    use = [header[0]] + sidm
    # Load only relevant cell-line columns; float32 footprint is manageable.
    df_parts = []
    with zipfile.ZipFile(tpm_zip) as z, z.open(member) as stream:
        chunks = pd.read_csv(
            stream, usecols=use, chunksize=CHUNK_GENES, low_memory=False
        )
        for ch in chunks:
            # First 3 lines in source are model_name, data_source and gene_symbol.
            gene = ch.iloc[:, 0].astype("string")
            ch = ch[~gene.isin(["model_name", "data_source", "gene_symbol"])].copy()
            if ch.empty:
                continue
            genes = ch.iloc[:, 0].astype("string").str.strip()
            matrix = (
                ch[sidm]
                .apply(pd.to_numeric, errors="coerce")
                .to_numpy(dtype=np.float32)
            )
            matrix[(matrix < 0) | (~np.isfinite(matrix))] = np.nan
            # Invalid all-missing rows will be excluded by training-only variance.
            df_parts.append((genes.to_numpy(dtype=str), matrix))
    if not df_parts:
        raise ValueError("No quantitative gene-expression rows found")
    genes = np.concatenate([p[0] for p in df_parts])
    vals = np.vstack([p[1] for p in df_parts])
    vals = np.log2(vals + 1.0)
    del df_parts
    tr_col = np.array([i for i, m in enumerate(sidm) if m in training_ids])
    if len(tr_col) < 100:
        raise ValueError("Fewer than 100 train cell lines with RNA")
    train_v = vals[:, tr_col]
    nonmissing = np.isfinite(train_v).mean(axis=1)
    means = np.nanmean(train_v, axis=1)
    variances = np.nanvar(train_v, axis=1)
    # Avoid missing, constant or duplicate gene labels; biologically valid annotations only.
    good = np.flatnonzero(
        (nonmissing >= 0.95)
        & np.isfinite(variances)
        & (variances > 0.01)
        & np.isfinite(means)
        & (means >= 0.1)
        & pd.Series(genes)
        .str.match(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", na=False)
        .to_numpy()
    )
    ranks = good[np.argsort(-variances[good], kind="stable")]
    selected = []
    seen = set()
    for i in ranks:
        symbol = genes[i].upper()
        if symbol not in seen:
            selected.append(i)
            seen.add(symbol)
        if len(selected) == N_GENES:
            break
    if len(selected) < 20:
        raise ValueError(f"Insufficient valid expression features: {len(selected)}")
    # cell lines x selected genes
    x = vals[selected, :].T.copy()
    chosen = [genes[i] for i in selected]
    info = {
        "total_gene_rows": len(genes),
        "expression_models": len(sidm),
        "candidate_gene_rows": len(good),
        "selected_genes": len(chosen),
        "missing_selected_values": int(np.isnan(x).sum()),
    }
    return pd.DataFrame(x, index=sidm, columns=chosen), info


def baseline_cancer_mean(train, other, key="cancer_type"):
    global_mean = float(train["AUC"].mean())
    avg = train.groupby(key)["AUC"].mean()
    return other[key].map(avg).fillna(global_mean).to_numpy(float)


def fit_models(xtrain, ytrain, xval, seed):
    """Predefined hyperparameters; validation used only for model selection."""
    models = {
        "ElasticNet": make_pipeline(
            SimpleImputer(strategy="median"),
            StandardScaler(),
            ElasticNet(alpha=0.03, l1_ratio=0.5, max_iter=3500, random_state=seed),
        ),
        "RandomForest": make_pipeline(
            SimpleImputer(strategy="median"),
            RandomForestRegressor(
                n_estimators=RF_ESTIMATORS,
                min_samples_leaf=5,
                max_features=0.5,
                n_jobs=-1,
                random_state=seed,
            ),
        ),
    }
    if HAS_XGB:
        models["XGBoost"] = make_pipeline(
            SimpleImputer(strategy="median"),
            XGBRegressor(
                n_estimators=XGB_ESTIMATORS,
                max_depth=3,
                learning_rate=0.04,
                subsample=0.8,
                colsample_bytree=0.8,
                min_child_weight=8,
                reg_lambda=8,
                tree_method="hist",
                n_jobs=4,
                random_state=seed,
            ),
        )
    fitted = {}
    for name, mdl in models.items():
        mdl.fit(xtrain, ytrain)
        fitted[name] = (mdl, mdl.predict(xval))
    return fitted


def main():
    if not DATA.is_dir():
        raise FileNotFoundError(f"Data folder unavailable: {DATA}")
    log(
        "SETUP",
        f"Python={sys.version.split()[0]} | pandas={pd.__version__} | XGBoost={HAS_XGB}",
    )
    log(
        "SETUP",
        f"Random seed={SEED}; maximum drugs={MAX_DRUGS}; selected genes={N_GENES}",
    )
    log(
        "DESIGN",
        "Cell line/patient grouped 70/15/15 TRAIN/VALIDATION/TEST, fixed before outcome modeling.",
    )
    log(
        "DESIGN",
        "Top drugs selected by TRAIN sample count only; high-variance genes selected using TRAIN expression only.",
    )
    log(
        "DESIGN",
        "All test sets remain unused until validation model choice is complete.",
    )
    models = pd.read_csv(get_file("model_list_"), low_memory=False)
    required = ["model_id", "patient_id", "tissue", "cancer_type"]
    if any(c not in models for c in required):
        raise ValueError(
            f"Missing model columns: {set(required) - set(models.columns)}"
        )
    models["model_id"] = norm(models["model_id"])
    models = models.dropna(subset=["model_id"]).drop_duplicates("model_id").copy()
    drugs = pd.read_excel(get_file("GDSC2_fitted_dose_response_"))
    req = ["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC", "LN_IC50"]
    if any(c not in drugs for c in req):
        raise ValueError(f"Missing GDSC2 columns: {set(req) - set(drugs.columns)}")
    drugs["model_id"] = norm(drugs["SANGER_MODEL_ID"])
    drugs["drug_id"] = norm(drugs["DRUG_ID"])
    drugs["AUC"] = pd.to_numeric(drugs["AUC"], errors="coerce")
    drugs["LN_IC50"] = pd.to_numeric(drugs["LN_IC50"], errors="coerce")
    drugs = drugs.dropna(subset=["model_id", "drug_id", "AUC"]).copy()
    drugs = drugs[np.isfinite(drugs["AUC"])].copy()
    drugs = drugs.merge(
        models[["model_id", "patient_id", "tissue", "cancer_type"]],
        on="model_id",
        how="inner",
        validate="many_to_one",
        suffixes=("", "_passport"),
    )
    if drugs.duplicated(["model_id", "drug_id"]).any():
        raise ValueError(
            "Duplicate model-drug response pairs: inspect source before continuing"
        )
    log(
        "GDSC",
        f"{len(drugs):,} valid records | {drugs.model_id.nunique():,} models | {drugs.drug_id.nunique():,} compounds",
    )
    # only models with RNA are eligible; first inspect ZIP header to identify columns.
    zpath = get_file("rnaseq_merged_")
    with zipfile.ZipFile(zpath) as z:
        m = [n for n in z.namelist() if "rsem_tpm_" in n.lower() and n.endswith(".csv")]
        if len(m) != 1:
            raise ValueError(f"Cannot locate unique wide TPM member: {m}")
        with z.open(m[0]) as stream:
            hdr = pd.read_csv(stream, nrows=0).columns
    eligible = set(drugs.model_id) & {str(c) for c in hdr if str(c).startswith("SIDM")}
    model_subset = (
        models[models.model_id.isin(eligible)]
        .copy()
        .sort_values("model_id")
        .reset_index(drop=True)
    )
    splits = pick_group_splits(model_subset)
    split_lookup = splits.set_index("model_id")["split"]
    drugs = drugs[drugs.model_id.isin(eligible)].copy()
    drugs["split"] = drugs.model_id.map(split_lookup)
    splits.to_csv(OUT / "model_splits.csv", index=False)
    split_ct = splits.split.value_counts()
    log(
        "SPLIT",
        f"Eligible with RNA: {len(eligible)} cell lines; "
        + "; ".join(f"{s}={split_ct.get(s, 0)}" for s in ("train", "val", "test")),
    )
    log(
        "SPLIT",
        "Actual response rows by split: " + str(drugs.split.value_counts().to_dict()),
    )
    log(
        "SPLIT",
        "Source patient_id grouping applied where available; exact cell-model groups never cross partitions.",
    )
    # Train-only ranking avoids future information leakage.
    c = (
        drugs[drugs.split == "train"]
        .groupby("drug_id")
        .agg(n=("model_id", "nunique"), auc_std=("AUC", "std"))
    )
    candidates = c.loc[(c.n >= MIN_TRAIN_PER_DRUG) & (c.auc_std >= 0.015)].sort_values(
        ["n", "auc_std"], ascending=False
    )
    selected_ids = candidates.head(MAX_DRUGS).index.tolist()
    if not selected_ids:
        raise ValueError(
            "No eligible drugs. Inspect training counts and eligibility thresholds."
        )
    label_map = drugs.groupby("drug_id").DRUG_NAME.first().to_dict()
    log(
        "DRUGS",
        f"{len(candidates)} eligible compounds; selected first {len(selected_ids)} by training coverage",
    )
    log(
        "DRUGS",
        "Selected: "
        + ", ".join(
            f"{label_map.get(i, i)} ({i}, train n={int(c.loc[i, 'n'])})"
            for i in selected_ids
        ),
    )
    dr = drugs[drugs.drug_id.isin(selected_ids)].copy()
    train_ids = set(splits.loc[splits.split == "train", "model_id"])
    expression, meta = read_expression(zpath, eligible, train_ids)
    log("RNA", str(meta))
    log(
        "RNA",
        "Scale: log2(TPM+1), not z-scores; ElasticNet standardisation uses training rows only.",
    )
    pd.DataFrame({"gene": expression.columns}).to_csv(
        OUT / "selected_genes.csv", index=False
    )
    # Store model-ready cell-level matrix; no raw source files modified.
    expression.rename_axis("model_id").reset_index().to_csv(
        OUT / "expression_selected_genes.csv.gz", index=False, compression="gzip"
    )
    all_validation = []
    all_test = []
    test_preds = []
    drug_coverage = []
    for idx, did in enumerate(selected_ids, 1):
        piece = dr[dr.drug_id == did].copy().sort_values("model_id")
        piece = piece[piece.model_id.isin(expression.index)].copy()
        tr = piece[piece.split == "train"].copy()
        va = piece[piece.split == "val"].copy()
        te = piece[piece.split == "test"].copy()
        drug_name = str(label_map.get(did, did))
        drug_coverage.append(
            {
                "drug_id": did,
                "drug_name": drug_name,
                "train": len(tr),
                "val": len(va),
                "test": len(te),
                "train_auc_std": float(tr.AUC.std()),
                "train_auc_median": float(tr.AUC.median()),
            }
        )
        if min(len(tr), len(va), len(te)) < 15:
            issue(
                "SKIP",
                f"{drug_name}: insufficient held-out data train/val/test={len(tr)}/{len(va)}/{len(te)}",
            )
            continue
        xtr, xva, xte = [
            expression.loc[rows.model_id].to_numpy(dtype=np.float32)
            for rows in (tr, va, te)
        ]
        ytr, yva, yte = [rows.AUC.to_numpy(dtype=float) for rows in (tr, va, te)]
        fitted = {}
        # Compare against honest training-only baselines (essential control).
        val_preds = {
            "TrainMean": np.full(len(va), ytr.mean()),
            "CancerMean": baseline_cancer_mean(tr, va),
        }
        test_fns = {
            "TrainMean": lambda: np.full(len(te), ytr.mean()),
            "CancerMean": lambda: baseline_cancer_mean(tr, te),
        }
        for name, (mdl, pred) in fit_models(xtr, ytr, xva, SEED + idx).items():
            fitted[name] = mdl
            val_preds[name] = pred
            test_fns[name] = lambda model=mdl: model.predict(xte)
        val_rows = []
        for name, pred in val_preds.items():
            row = {
                "drug_id": did,
                "drug_name": drug_name,
                "model": name,
                **metrics(yva, pred),
            }
            val_rows.append(row)
            all_validation.append(row)
        # Selection uses validation MAE; never use test scores to pick a winner.
        val_sorted = sorted(val_rows, key=lambda r: (r["MAE"], r["model"]))
        winner = val_sorted[0]["model"]
        for name in ("TrainMean", "CancerMean", winner):
            pred = np.asarray(test_fns[name](), dtype=float)
            test_row = {
                "drug_id": did,
                "drug_name": drug_name,
                "model": name,
                "validation_selected": name == winner,
                **metrics(yte, pred),
            }
            all_test.append(test_row)
            if name == winner:
                test_preds.extend(
                    {
                        "model_id": mid,
                        "drug_id": did,
                        "drug_name": drug_name,
                        "cancer_type": ct,
                        "observed_AUC": float(obs),
                        "predicted_AUC": float(p),
                    }
                    for mid, ct, obs, p in zip(te.model_id, te.cancer_type, yte, pred)
                )
        log(
            "MODEL",
            f"{idx:02d}/{len(selected_ids)} {drug_name} | n={len(tr)}/{len(va)}/{len(te)} "
            f"| selected={winner} | val MAE={val_sorted[0]['MAE']:.4f} "
            f"| test MAE={next(t['MAE'] for t in all_test[::-1] if t['drug_id'] == did and t['model'] == winner):.4f}",
        )
    pd.DataFrame(drug_coverage).to_csv(OUT / "drug_coverage.csv", index=False)
    pd.DataFrame(all_validation).to_csv(
        OUT / "validation_model_comparison.csv", index=False
    )
    pd.DataFrame(all_test).to_csv(OUT / "test_model_comparison.csv", index=False)
    pd.DataFrame(test_preds).to_csv(OUT / "heldout_test_predictions.csv", index=False)
    if not all_test:
        raise RuntimeError("No models could be evaluated. See warnings.")
    test_df = pd.DataFrame(all_test)
    sel = test_df[test_df.validation_selected].copy()
    base = test_df[test_df.model == "CancerMean"].set_index("drug_id")
    sel = format_selected_results(sel, base)
    sel.sort_values("MAE_improvement", ascending=False).to_csv(
        OUT / "selected_model_vs_baseline.csv", index=False
    )
    log(
        "RESULT",
        f"Modelled drugs: {len(sel)}; validation-selected model beats cancer-mean on held-out TEST: "
        f"{int(sel.improved.sum())}/{len(sel)}",
    )
    log(
        "RESULT",
        f"Median TEST MAE selected={sel.MAE.median():.4f}; cancer-mean={sel.baseline_MAE.median():.4f}; "
        f"median improvement={sel.MAE_improvement.median():+.4f}",
    )
    log("RESULT", "Selected model types: " + str(sel.model.value_counts().to_dict()))
    log(
        "RESULT",
        f"Median TEST R2={sel.R2.median():.4f}, Spearman={sel.Spearman.median():.4f}",
    )
    for r in sel.sort_values(
        "MAE_improvement", ascending=False
    ).itertuples():
        log(
            "PER-DRUG",
            f"{r.drug_name} | model={r.model} | test n={r.n} | MAE={r.MAE:.4f} | "
            f"cancer baseline MAE={r.baseline_MAE:.4f} | delta={r.MAE_improvement:+.4f} "
            f"| R2={r.R2:.3f} | Spearman={r.Spearman:.3f}",
        )
    log(
        "LIMITATION",
        "Pilot assesses 15 high-coverage drugs, not all 295; GNN, PRISM, mutation and proteomics belong to Phase 3/4.",
    )
    log(
        "LIMITATION",
        "An unseen-cell-line test does NOT prove prediction on unseen drugs, cancer tissues or patients.",
    )
    log(
        "LIMITATION",
        "TPM harmonised from Sanger/Broad sources; batch/source effects and assay comparability need later sensitivity analysis.",
    )
    log(
        "LIMITATION",
        "Testing selection across 15 drugs is exploratory; treat aggregate patterns as descriptive, not confirmatory.",
    )
    log("COMPLETION", f"Output directory: {OUT}")
    log(
        "COMPLETION",
        f"Runtime: {(time.time() - START) / 60:.1f} minutes; warnings: {len(ISSUES)}",
    )


def finish():
    text = "\n".join(
        [
            "=" * 82,
            "PROJECT 5 | RNA MODEL RESULTS",
            "=== BEGIN PHASE 2 RESULTS ===",
            *REPORT,
            f"[STATUS] Warnings/errors: {len(ISSUES)}",
            "=== END PHASE 2 RESULTS ===",
            "=" * 82,
        ]
    )
    (OUT / "phase2_report.txt").write_text(text, encoding="utf-8")
    print("\n" + text, flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        import traceback

        issue("FATAL", f"{type(exc).__name__}: {exc}")
        log("TRACEBACK", traceback.format_exc()[-3000:])
        raise
    finally:
        finish()
