"""Evaluate modality performance with fold-specific feature selection.

Imputation, scaling and selection are learned inside each training fold, with identical held-out observations used for modality comparisons.

Run from the repository root with ``python 06_nested_validation.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# Fold-specific multi-omics validation
# ----------------------------------------------------------------------------

from pathlib import Path
import importlib.util, zipfile, warnings, traceback, time
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score
from xgboost import XGBRegressor

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
P2 = ROOT / "Results" / "Phase2"
OUT = ROOT / "Results" / "Phase5B"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261028
REPEATS = 3
FOLDS = 4
COUNTS = {"RNA": 256, "PRO": 256, "CNV": 128, "MUT": 64}
LOG = []
START = time.time()


def log(k, v):
    LOG.append(f"[{k}] {v}")


def norm(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def datafile(prefix):
    candidates = list(DATA.glob(prefix + "*"))
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected one Data/{prefix}*, found {[p.name for p in candidates]}"
        )
    return candidates[0]


def phase3():
    file = ROOT / "04_multiomics_explainability.py"
    if not file.exists():
        raise FileNotFoundError(f"Missing earlier project script {file}")
    spec = importlib.util.spec_from_file_location("p5_phase3a_helpers", file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_all_rna(models):
    """Use original wide RSEM TPM table, not the previously selected 256-gene matrix."""
    path = datafile("rnaseq_merged_")
    with zipfile.ZipFile(path) as z:
        members = [
            n
            for n in z.namelist()
            if n.lower().endswith(".csv") and "rsem_tpm_" in n.lower()
        ]
        if len(members) != 1:
            raise ValueError(f"Unexpected RNA zip CSV members: {members}")
        with z.open(members[0]) as h:
            header = pd.read_csv(h, nrows=0).columns.tolist()
        cols = [c for c in header if c.startswith("SIDM") and c in models]
        if len(cols) < 300:
            raise ValueError(f"Only {len(cols)} RNA-matched training cells")
        gene_columns = [header[0]]
        parts = []
        labels = []
        with z.open(members[0]) as h:
            for chunk in pd.read_csv(
                h, usecols=gene_columns + cols, chunksize=1000, low_memory=False
            ):
                sym = chunk.iloc[:, 0].astype("string").str.strip()
                mask = (
                    ~sym.isin(["model_name", "data_source", "gene_symbol"])
                    & sym.notna()
                )
                if not mask.any():
                    continue
                chunk = chunk.loc[mask]
                labels.extend(sym.loc[mask].astype(str).tolist())
                parts.append(
                    chunk[cols]
                    .apply(pd.to_numeric, errors="coerce")
                    .to_numpy(dtype=np.float32)
                )
    matrix = np.vstack(parts)
    matrix[~np.isfinite(matrix) | (matrix < 0)] = np.nan
    np.log2(matrix + 1, out=matrix)
    # Exclude ambiguous duplicate symbols globally based only on annotations, not outcomes.
    symbols = pd.Index([s.upper() for s in labels])
    valid = symbols.str.match(r"^[A-Z0-9][A-Z0-9_.-]*$") & ~symbols.duplicated(
        keep=False
    )
    matrix = matrix[valid, :]
    labels = symbols[valid]
    log(
        "RNA",
        f"Full RNA loaded: {len(labels):,} uniquely named genes and {len(cols)} original TRAIN models; fold-wise gene selection",
    )
    return pd.DataFrame(matrix.T, index=cols, columns=labels)


def rank_continuous(
    mat, indices, k, minimum_coverage=0.70, minimum_variance=1e-8, minimum_mean=None
):
    x = mat.iloc[indices]
    coverage = x.notna().mean(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        variances = x.var(axis=0).replace([np.inf, -np.inf], np.nan).fillna(0)
        means = x.mean(axis=0).replace([np.inf, -np.inf], np.nan)
    valid = (coverage >= minimum_coverage) & (variances > minimum_variance)
    if minimum_mean is not None:
        valid &= means >= minimum_mean
    return variances[valid].nlargest(k).index.tolist()


def feature_sets(mats, train_idx):
    selected = {
        "RNA": rank_continuous(mats["RNA"], train_idx, 256, 0.95, 0.01, 0.1),
        "PRO": rank_continuous(mats["PRO"], train_idx, 256, 0.70),
        "CNV": rank_continuous(mats["CNV"], train_idx, 128, 0.70),
    }
    mut = mats["MUT"].iloc[train_idx].fillna(0)
    counts = mut.sum(axis=0)
    selected["MUT"] = (
        counts[(counts >= 3) & (counts <= len(train_idx) - 3)]
        .nlargest(64)
        .index.tolist()
    )
    if any(len(cols) < min(20, COUNTS[k]) for k, cols in selected.items()):
        raise RuntimeError(
            f"Insufficient predictors inside fold: { {k: len(v) for k, v in selected.items()} }"
        )
    return selected


def fit_predict(mats, sel, ti, vi, y, kind, seed):
    blocks = {
        "RNA": ["RNA"],
        "Genomic": ["RNA", "MUT", "CNV"],
        "AllFour": ["RNA", "PRO", "MUT", "CNV"],
    }
    train = np.hstack(
        [mats[b].iloc[ti][sel[b]].to_numpy(dtype=np.float32) for b in blocks[kind]]
    )
    val = np.hstack(
        [mats[b].iloc[vi][sel[b]].to_numpy(dtype=np.float32) for b in blocks[kind]]
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Skipping features without any observed values"
        )
        imp = SimpleImputer(strategy="median", keep_empty_features=True)
        train = imp.fit_transform(train)
        val = imp.transform(val)
    reg = XGBRegressor(
        n_estimators=240,
        max_depth=3,
        learning_rate=0.045,
        subsample=0.85,
        colsample_bytree=0.85,
        min_child_weight=8,
        reg_lambda=8,
        objective="reg:squarederror",
        tree_method="hist",
        n_jobs=4,
        random_state=seed,
        verbosity=0,
    )
    reg.fit(train, y[ti])
    return reg.predict(val)


def main():
    log(
        "DESIGN",
        "Original TRAIN only; fold-wise RNA/protein/CNV/mutation selection; patient-group-aware 3x4 CV.",
    )
    log(
        "DESIGN",
        "Identical evaluation models per modality, fold-fitted imputation, fixed hyperparameters; no VAL/TEST responses accessed.",
    )
    splits = pd.read_csv(P2 / "model_splits.csv", dtype={"model_id": "string"})
    splits.model_id = norm(splits.model_id)
    if splits.model_id.duplicated().any():
        raise ValueError("Duplicate model_splits model IDs")
    if "group_key" not in splits.columns:
        raise ValueError(
            "Original group_key required to prevent derived-patient leakage"
        )
    eligible = set(splits.loc[splits.split.eq("train"), "model_id"])
    p3 = phase3()
    # Raw omics data, without using previously selected feature lists.
    rna = read_all_rna(eligible)
    prot = p3.read_proteome(eligible)
    cnv = p3.read_cnv(eligible)
    mut = p3.read_mutation(eligible)
    av = pd.read_csv(datafile("model_dataset_availability_"), low_memory=False)
    av.model_id = norm(av.model_id)
    if av.model_id.duplicated().any():
        raise ValueError("Duplicated model availability IDs")
    av = av.set_index("model_id")
    covcols = [
        c
        for c in ("Mutation Sanger Cell Lines (WES)", "Mutation Broad Cell Lines (WES)")
        if c in av
    ]
    if not covcols:
        raise ValueError("No mutation coverage metadata available")
    covered = pd.Series(False, index=av.index)
    for c in covcols:
        covered |= av[c].astype("string").str.casefold().eq("true").fillna(False)
    matched = sorted(
        eligible
        & set(rna.index)
        & set(prot.index)
        & set(cnv.index)
        & set(covered[covered].index)
    )
    if len(matched) < 400:
        raise ValueError(f"Unexpectedly low matched count {len(matched)}")
    mats = {
        "RNA": rna.reindex(matched).astype("float32"),
        "PRO": prot.reindex(matched).astype("float32"),
        "CNV": cnv.reindex(matched).astype("float32"),
        "MUT": mut.reindex(matched).fillna(0).astype("float32"),
    }
    log(
        "COHORT",
        f"Matched original TRAIN: {len(matched)} cell lines; features before fold selection "
        + str({k: v.shape[1] for k, v in mats.items()}),
    )
    meta = splits.set_index("model_id").loc[matched]
    cancer = meta.cancer_type.astype("string")
    grouping = meta.group_key.astype("string")
    if grouping.isna().any():
        raise ValueError("Missing patient group identifiers")
    raw = pd.read_excel(
        datafile("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    raw["model_id"] = norm(raw.SANGER_MODEL_ID)
    raw["drug_id"] = norm(raw.DRUG_ID)
    raw["AUC"] = pd.to_numeric(raw.AUC, errors="coerce")
    raw = raw[raw.model_id.isin(matched) & raw.AUC.between(0, 1)]
    if raw.duplicated(["model_id", "drug_id"]).any():
        raise ValueError("Duplicated drug-model measurements")
    d = pd.read_csv(P2 / "validation_model_comparison.csv", dtype={"drug_id": "string"})
    dids = norm(d.drug_id).dropna().drop_duplicates().tolist()
    if len(dids) != 15:
        raise ValueError(f"Expected 15 selected drugs, got {len(dids)}")
    allrows = []
    fold_features = []
    oofrows = []
    for ix, did in enumerate(dids, 1):
        sub = raw[raw.drug_id.eq(did)].set_index("model_id")
        ids = sorted(set(matched) & set(sub.index))
        if len(ids) < 200:
            raise ValueError(f"{did} has only {len(ids)} eligible training lines")
        index = pd.Index(matched).get_indexer(ids)
        y = sub.loc[ids, "AUC"].to_numpy(dtype=float)
        target_can = cancer.reindex(ids).to_numpy()
        target_group = grouping.reindex(ids).astype(str).to_numpy()
        unique_groups = np.unique(target_group)
        drug_name = str(sub.DRUG_NAME.iloc[0])
        drug_rows = []
        for rep in range(REPEATS):
            # Deterministic random group label assignment, fed to GroupKFold.
            rng = np.random.default_rng(SEED + rep)
            shuffle = rng.permutation(len(unique_groups))
            mapped = dict(zip(unique_groups, shuffle))
            randomized = np.array([mapped[g] for g in target_group])
            splitter = GroupKFold(n_splits=FOLDS)
            for fold, (ti, vi) in enumerate(splitter.split(ids, groups=randomized), 1):
                if set(target_group[ti]) & set(target_group[vi]):
                    raise ValueError("Patient-group leakage")
                # Restrict all omics to drug-specific TRAIN IDs; full feature ranking occurs here.
                x = {k: v.iloc[index].reset_index(drop=True) for k, v in mats.items()}
                chosen = feature_sets(x, ti)
                for key, features in chosen.items():
                    fold_features.extend(
                        {
                            "drug_id": did,
                            "repeat": rep + 1,
                            "fold": fold,
                            "modality": key,
                            "feature": f,
                        }
                        for f in features
                    )
                mean = float(y[ti].mean())
                cmeans = (
                    pd.DataFrame({"cancer": target_can[ti], "y": y[ti]})
                    .groupby("cancer")
                    .y.mean()
                )
                baseline = (
                    pd.Series(target_can[vi]).map(cmeans).fillna(mean).to_numpy(float)
                )
                models = {"TrainMean": np.repeat(mean, len(vi)), "CancerMean": baseline}
                for kind in ["RNA", "Genomic", "AllFour"]:
                    models[kind] = fit_predict(
                        x, chosen, ti, vi, y, kind, SEED + rep * 100 + fold
                    )
                for kind, pred in models.items():
                    mae = float(mean_absolute_error(y[vi], pred))
                    r2 = (
                        float(r2_score(y[vi], pred))
                        if np.std(y[vi]) > 1e-10
                        else np.nan
                    )
                    drug_rows.append(
                        dict(
                            drug_id=did,
                            drug_name=drug_name,
                            repeat=rep + 1,
                            fold=fold,
                            modality=kind,
                            n_train=len(ti),
                            n_valid=len(vi),
                            MAE=mae,
                            R2=r2,
                        )
                    )
                    for index_val, y_true, y_hat in zip(vi, y[vi], pred):
                        oofrows.append(
                            dict(
                                drug_id=did,
                                drug_name=drug_name,
                                repeat=rep + 1,
                                fold=fold,
                                model_id=ids[index_val],
                                modality=kind,
                                observed_auc=float(y_true),
                                predicted_auc=float(y_hat),
                            )
                        )
        allrows.extend(drug_rows)
        log(
            "DRUG",
            f"{ix:02d}/15 {drug_name} [{did}] n={len(ids)}; patient-grouped folds={REPEATS * FOLDS}",
        )
    scores = pd.DataFrame(allrows)
    scores.to_csv(OUT / "foldwise_cv_model_scores.csv", index=False)
    pd.DataFrame(fold_features).to_csv(
        OUT / "foldwise_selected_feature_lists.csv", index=False
    )
    preds = pd.DataFrame(oofrows)
    preds.to_csv(OUT / "foldwise_oof_predictions.csv", index=False)
    paired = scores.pivot(
        index=["drug_id", "drug_name", "repeat", "fold"],
        columns="modality",
        values="MAE",
    ).reset_index()
    paired["gain_all_vs_rna"] = paired.RNA - paired.AllFour
    paired["gain_all_vs_cancer"] = paired.CancerMean - paired.AllFour
    paired.to_csv(OUT / "foldwise_paired_performance.csv", index=False)
    drug = paired.groupby(["drug_id", "drug_name"], as_index=False).agg(
        n_folds=("fold", "size"),
        gain_vs_rna=("gain_all_vs_rna", "mean"),
        gain_vs_cancer=("gain_all_vs_cancer", "mean"),
        better_folds=("gain_all_vs_rna", lambda x: int((x > 0).sum())),
    )
    drug = drug.sort_values("gain_vs_rna", ascending=False)
    drug.to_csv(OUT / "foldwise_drug_summary.csv", index=False)
    summary = scores.groupby("modality", as_index=False).agg(
        mean_MAE=("MAE", "mean"),
        median_MAE=("MAE", "median"),
        median_R2=("R2", "median"),
    )
    summary.to_csv(OUT / "foldwise_modality_summary.csv", index=False)
    # One prediction per cell line per repeat; use repeat-level OOF MAE as descriptive stability.
    repeat_mae = (
        preds.groupby(["drug_id", "drug_name", "repeat", "modality"])
        .apply(
            lambda g: mean_absolute_error(g.observed_auc, g.predicted_auc),
            include_groups=False,
        )
        .rename("OOF_MAE")
        .reset_index()
    )
    repeat_mae.to_csv(OUT / "repeat_level_oof_MAE.csv", index=False)
    # Drug-level bootstrap: descriptive uncertainty for 15 prespecified compounds, NOT population inference.
    gains = drug.gain_vs_rna.to_numpy()
    rng = np.random.default_rng(SEED)
    draws = np.mean(rng.choice(gains, size=(10000, len(gains)), replace=True), axis=1)
    ci = np.quantile(draws, [0.025, 0.975])
    log(
        "SUMMARY",
        "\n"
        + summary.sort_values("mean_MAE").to_string(
            index=False, float_format=lambda a: f"{a:.5f}"
        ),
    )
    log(
        "PAIRED",
        f"AllFour lower mean MAE vs RNA for {(drug.gain_vs_rna > 0).sum()}/15 drugs; mean paired drug gain={gains.mean():+.6f}; descriptive drug bootstrap 95% interval=[{ci[0]:+.6f},{ci[1]:+.6f}]",
    )
    log(
        "DRUG RESULTS",
        "\n"
        + drug[
            ["drug_name", "n_folds", "gain_vs_rna", "gain_vs_cancer", "better_folds"]
        ].to_string(index=False, float_format=lambda a: f"{a:.6f}"),
    )
    log(
        "REPEAT OOF",
        "\n"
        + repeat_mae[
            repeat_mae.drug_name.isin(["Nutlin-3a (-)", "Camptothecin", "PD0325901"])
        ]
        .pivot_table(
            index=["drug_name", "repeat"], columns="modality", values="OOF_MAE"
        )
        .reset_index()
        .to_string(index=False, float_format=lambda a: f"{a:.5f}"),
    )
    log(
        "CAUTION",
        "Selection of both features and imputation is inside training folds; patient-group independence checked. Fixed settings, preselected drugs, original TRAIN/assay-coverage cohort; no new independent TEST.",
    )
    log(
        "CAUTION",
        "Drug bootstrap interval describes variability among the 15 chosen drugs, not sampling precision over all compounds or all clinical cancers. Repeated folds must not be treated as independent replicates.",
    )
    log(
        "STATUS",
        f"SUCCESS; elapsed={(time.time() - START) / 60:.1f} min; outputs={OUT}",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log("STATUS", f"FAILED: {type(exc).__name__}: {exc}")
        log("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 94 + "\n=== BEGIN PROJECT 5 PHASE 5B RESULTS ===")
        print("\n".join(LOG))
        print("=== END PROJECT 5 PHASE 5B RESULTS ===\n" + "=" * 94)
