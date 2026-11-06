"""Compare RNA and multi-omics prediction and calculate model explanations.

SHAP values are treated as predictive attributions rather than evidence of causal molecular mechanisms.

Run from the repository root with ``python 04_multiomics_explainability.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# Multi-omics predictive comparison
# ----------------------------------------------------------------------------

from pathlib import Path
import zipfile, traceback, time, sys, warnings
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
P2 = ROOT / "Results" / "Phase2"
OUT = ROOT / "Results" / "Phase3A"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261028
N_PROTEINS = 256
N_CNV = 128
N_MUTATIONS = 64
MIN_TRAIN = 120
MIN_VAL = 25
REPORT = []
T0 = time.time()


def log(k, msg):
    REPORT.append(f"[{k}] {msg}")


def norm(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def file(prefix):
    a = list(DATA.glob(prefix + "*"))
    if len(a) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {prefix}*, found {[p.name for p in a]}"
        )
    return a[0]


def zipcsv(path, fragment, **kw):
    with zipfile.ZipFile(path) as z:
        a = [
            x
            for x in z.namelist()
            if fragment.lower() in Path(x).name.lower() and x.lower().endswith(".csv")
        ]
        if len(a) != 1:
            raise ValueError(f"Cannot uniquely choose CSV containing {fragment}: {a}")
        with z.open(a[0]) as stream:
            return pd.read_csv(stream, **kw)


def fmt(df, cols, n=20):
    if df.empty:
        return "No eligible rows"
    return df[cols].head(n).to_string(index=False, float_format=lambda x: f"{x:.5g}")


def choose_features(mat, training_ids, k, kind, min_coverage=0.65):
    """Screen using TRAIN molecular data only; missingness and variance measured on TRAIN."""
    tr = mat.reindex(training_ids)
    cov = tr.notna().mean(axis=0)
    var = tr.var(axis=0, skipna=True).replace([np.inf, -np.inf], np.nan).fillna(0)
    eligible = (cov >= min_coverage) & (var > 1e-10)
    chosen = var[eligible].sort_values(ascending=False).head(k).index.tolist()
    log(
        "FEATURES",
        f"{kind}: candidates={mat.shape[1]:,}; train coverage>={min_coverage:.0%} and variable={int(eligible.sum()):,}; chosen={len(chosen)}",
    )
    return chosen


def read_proteome(models):
    path = file("Proteomics_")
    with zipfile.ZipFile(path) as z:
        names = [
            n
            for n in z.namelist()
            if Path(n).name.lower().startswith("proteomics_all_") and n.endswith(".csv")
        ]
        if len(names) != 1:
            raise ValueError(f"Proteomics long table ambiguous: {names}")
        with z.open(names[0]) as h:
            hdr = pd.read_csv(h, nrows=0).columns.tolist()
        log("SCHEMA", f"Proteomics source columns: {hdr}")
        required = {"model_id", "symbol", "protein_intensity"}
        if not required.issubset(hdr):
            raise ValueError(
                f"Proteomics missing expected columns: {required - set(hdr)}"
            )
        chunks = []
        with z.open(names[0]) as h:
            for ch in pd.read_csv(
                h, usecols=list(required), chunksize=300000, low_memory=False
            ):
                ch.model_id = norm(ch.model_id)
                ch = ch[ch.model_id.isin(models) & ch.symbol.notna()].copy()
                ch["protein_intensity"] = pd.to_numeric(
                    ch.protein_intensity, errors="coerce"
                )
                ch = ch[np.isfinite(ch.protein_intensity)]
                ch.symbol = ch.symbol.astype(str).str.strip()
                ch = ch[ch.symbol.ne("")]
                if not ch.empty:
                    chunks.append(ch)
    if not chunks:
        raise ValueError("No proteomics measurements aligned to Phase2 models")
    long = pd.concat(chunks, ignore_index=True)
    del chunks
    # Mean of same-symbol quantitative records per model is explicitly documented,
    # not silently treating multiple assays/protein groups as independent predictors.
    dups = int(long.duplicated(["model_id", "symbol"]).sum())
    log(
        "PROTEOME",
        f"Rows={len(long):,}; unique models={long.model_id.nunique()}; duplicate model-symbol records={dups:,}; mean aggregation for duplicates",
    )
    wide = (
        long.groupby(["model_id", "symbol"], observed=True)
        .protein_intensity.mean()
        .unstack("symbol")
        .astype("float32")
    )
    return wide


def read_cnv(models):
    df = zipcsv(
        file("cnv_summary_"),
        "cnv_summary_",
        usecols=["model_id", "symbol", "total_copy_number", "source", "data_type"],
        low_memory=False,
    )
    df.model_id = norm(df.model_id)
    df = df[df.model_id.isin(models) & df.symbol.notna()].copy()
    df["total_copy_number"] = pd.to_numeric(df.total_copy_number, errors="coerce")
    df = df[np.isfinite(df.total_copy_number) & df.total_copy_number.between(0, 100)]
    df.symbol = df.symbol.astype(str).str.strip()
    log(
        "CNV",
        f"Raw matching records={len(df):,}; assay sources={df.source.value_counts(dropna=False).head(6).to_dict()}; data types={df.data_type.value_counts(dropna=False).head(6).to_dict()}",
    )
    # Do not silently combine records from different assays. Keep only genes with
    # a single distinct quantitative record per cell line; other combinations missing.
    cnt = df.groupby(["model_id", "symbol"]).size().rename("n").reset_index()
    df = df.merge(cnt, on=["model_id", "symbol"], validate="many_to_one")
    duplicate_n = int((df.n > 1).sum())
    df = df[df.n == 1]
    log(
        "CNV",
        f"Excluded {duplicate_n:,} rows from duplicated model-gene/assay combinations; no assay-scale pooling",
    )
    return df.pivot(
        index="model_id", columns="symbol", values="total_copy_number"
    ).astype("float32")


def read_mutation(models):
    df = zipcsv(
        file("mutations_summary_"),
        "mutations_summary_",
        usecols=["model_id", "gene_symbol"],
        low_memory=False,
    )
    df.model_id = norm(df.model_id)
    df = df[df.model_id.isin(models) & df.gene_symbol.notna()].copy()
    df.gene_symbol = df.gene_symbol.astype(str).str.strip()
    df = df[df.gene_symbol.ne("")].drop_duplicates(["model_id", "gene_symbol"])
    log(
        "MUTATION",
        f"Annotated mutation events={len(df):,}; observed models={df.model_id.nunique()}; this is driver mutation SUMMARY, not whole-exome variant calls",
    )
    if df.empty:
        raise ValueError("No matched mutation records")
    x = pd.crosstab(df.model_id, df.gene_symbol).clip(upper=1).astype("float32")
    # Zero = no event in summary, NEVER necessarily reference/wild type.
    return x


def trial_model(xtr, ytr, xva, yva, tag):
    # Fit missing-value and scaling transformations on TRAIN only.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Skipping features without any observed values"
        )
        fill = SimpleImputer(strategy="median", keep_empty_features=True)
        xtr = fill.fit_transform(xtr)
        xva = fill.transform(xva)
    scaler = StandardScaler()
    xtr = scaler.fit_transform(xtr).astype("float32")
    xva = scaler.transform(xva).astype("float32")
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
        random_state=SEED,
        verbosity=0,
    )
    reg.fit(xtr, ytr)
    yp = reg.predict(xva)
    return float(mean_absolute_error(yva, yp)), float(r2_score(yva, yp)), yp, reg


def main():
    log(
        "DESIGN",
        "Matched-cohort comparison across 15 GDSC2 drug IDs; all models use identical TRAIN and VAL samples per compound.",
    )
    log(
        "DESIGN",
        "Original grouped split preserved; TEST untouched (previously inspected in Phase2). Validation reuse is exploratory, NOT an independent confirmatory estimate.",
    )
    log(
        "DESIGN",
        "XGBoost fixed parameters for every modality; protein/CNV/mutation features ranked on TRAIN only. Median imputation and scaling fitted on TRAIN only.",
    )
    log(
        "DESIGN",
        "Primary metric MAE on GDSC2 AUC; lower MAE better. No cross-drug AUC comparisons or causal/clinical interpretations.",
    )
    splits = pd.read_csv(P2 / "model_splits.csv")
    splits.model_id = norm(splits.model_id)
    if splits.model_id.duplicated().any():
        raise ValueError("Duplicate split IDs")
    if not {"train", "val", "test"}.issuperset(set(splits.split)):
        raise ValueError("Unrecognised split")
    rna = pd.read_csv(
        P2 / "expression_selected_genes.csv.gz", compression="gzip", low_memory=False
    )
    rna.model_id = norm(rna.model_id)
    if rna.model_id.duplicated().any():
        raise ValueError("Duplicate RNA models")
    rna = rna.set_index("model_id").apply(pd.to_numeric, errors="coerce")
    if rna.shape[1] != 256:
        log(
            "CAUTION",
            f"Found {rna.shape[1]} Phase2 RNA predictors instead of expected 256",
        )
    ids = set(rna.index) & set(splits.model_id)
    pro = read_proteome(ids)
    cnv = read_cnv(ids)
    mut = read_mutation(ids)
    # This study deliberately limits to a complete assay-coverage cohort, to avoid
    # comparing models on different cell lines. Missing individual feature values
    # are handled with train-fitted imputation, not sample exclusion.
    available = pd.read_csv(file("model_dataset_availability_"), low_memory=False)
    available.model_id = norm(available.model_id)
    av = available.drop_duplicates("model_id").set_index("model_id")
    mut_cov_cols = [
        c
        for c in av
        if c in ("Mutation Sanger Cell Lines (WES)", "Mutation Broad Cell Lines (WES)")
    ]
    if not mut_cov_cols:
        raise ValueError("No mutation assay coverage columns in availability file")
    covered = pd.Series(False, index=av.index)
    for c in mut_cov_cols:
        covered |= av[c].astype("string").str.lower().eq("true").fillna(False)
    matched = sorted(
        ids & set(pro.index) & set(cnv.index) & set(covered[covered].index)
    )
    log(
        "COHORT",
        f"Original split models={len(ids)}; proteome models={len(set(pro.index) & ids)}; CNV models={len(set(cnv.index) & ids)}; mutation-assay-covered={len(set(covered[covered].index) & ids)}; matched={len(matched)}",
    )
    if len(matched) < 150:
        raise ValueError(
            f"Unexpectedly small matched cohort ({len(matched)}), inspect coverage metadata and source assays"
        )
    s = splits.set_index("model_id").loc[matched, "split"]
    trids = s[s == "train"].index.tolist()
    vaids = s[s == "val"].index.tolist()
    log(
        "SPLIT",
        f"matched TRAIN={len(trids)}, VAL={len(vaids)}, TEST={sum(s == 'test')} (TEST not accessed in modelling)",
    )
    # Original 15 selected compounds, as defined by Phase2 training-coverage selection.
    selected = (
        pd.read_csv(P2 / "drug_coverage.csv")
        if (P2 / "drug_coverage.csv").exists()
        else None
    )
    drug = pd.read_excel(
        file("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    drug["model_id"] = norm(drug.SANGER_MODEL_ID)
    drug["drug_id"] = norm(drug.DRUG_ID)
    drug.AUC = pd.to_numeric(drug.AUC, errors="coerce")
    drug = drug[drug.model_id.isin(matched) & drug.AUC.between(0, 1)].copy()
    if drug.duplicated(["model_id", "drug_id"]).any():
        raise ValueError("Duplicate drug-model observations")
    # Fixed TRAIN-ranked drug selection as originally performed, based only on original TRAIN.
    if selected is not None:
        log("SCHEMA", f"drug_coverage.csv columns={selected.columns.tolist()}")
    # Use Phase2 validation model comparison identifiers, containing original 15, rather than selecting on Phase3 outcomes.
    original = pd.read_csv(P2 / "validation_model_comparison.csv")
    log("SCHEMA", f"Phase2 validation table columns={original.columns.tolist()}")
    idcol = next((c for c in ["drug_id", "DRUG_ID"] if c in original), None)
    if idcol is None:
        raise ValueError(
            "Phase2 validation comparison lacks drug_id; cannot recover prespecified 15 drugs"
        )
    dids = norm(original[idcol]).dropna().unique().tolist()
    if len(dids) != 15:
        raise ValueError(f"Expected 15 unique original drugs, found {len(dids)}")
    # Training-selected feature sets independent of compound response.
    pc = choose_features(pro, trids, N_PROTEINS, "protein intensity", 0.70)
    cc = choose_features(cnv, trids, N_CNV, "CNV copy number", 0.70)
    # mutation summary: only event presence >2 and <= 50% among TRAIN; no-record treated as uncertainty.
    mut = mut.reindex(index=matched, fill_value=0)
    mprev = mut.reindex(trids, fill_value=0).sum(axis=0)
    mc = (
        mprev[(mprev >= 3) & (mprev <= len(trids) - 3)]
        .sort_values(ascending=False)
        .head(N_MUTATIONS)
        .index.tolist()
    )
    log(
        "FEATURES",
        f"Mutation event/no-record predictors selected={len(mc)}; no-record DOES NOT imply wild type",
    )
    if not pc or not cc or not mc:
        raise ValueError(
            "A modality has zero eligible features; inspect logged source schema"
        )
    # CNV baseline magnitude varies between assay sources. Keep interpretative caution in report.
    features = {
        "RNA": rna.reindex(matched).astype("float32"),
        "Protein": pro.reindex(matched)[pc].astype("float32"),
        "RNA+Protein": pd.concat(
            [
                rna.reindex(matched).add_prefix("RNA:"),
                pro.reindex(matched)[pc].add_prefix("PRO:"),
            ],
            axis=1,
        ).astype("float32"),
        "RNA+Mutation+CNV": pd.concat(
            [
                rna.reindex(matched).add_prefix("RNA:"),
                mut.reindex(matched)[mc].add_prefix("MUT:"),
                cnv.reindex(matched)[cc].add_prefix("CNV:"),
            ],
            axis=1,
        ).astype("float32"),
        "AllFour": pd.concat(
            [
                rna.reindex(matched).add_prefix("RNA:"),
                pro.reindex(matched)[pc].add_prefix("PRO:"),
                mut.reindex(matched)[mc].add_prefix("MUT:"),
                cnv.reindex(matched)[cc].add_prefix("CNV:"),
            ],
            axis=1,
        ).astype("float32"),
    }
    selected_feats = pd.concat(
        [
            pd.DataFrame({"modality": kind, "feature": cols})
            for kind, cols in [("protein", pc), ("CNV", cc), ("mutation", mc)]
        ],
        ignore_index=True,
    )
    selected_feats.to_csv(OUT / "selected_multiomics_features.csv", index=False)
    meta = pd.read_csv(
        file("model_list_"), usecols=["model_id", "cancer_type"], low_memory=False
    )
    meta.model_id = norm(meta.model_id)
    cancer = meta.drop_duplicates("model_id").set_index("model_id").cancer_type
    records = []
    predictions = []
    for did in dids:
        ds = drug[drug.drug_id == did].set_index("model_id")
        tids = ds.index.intersection(trids)
        vids = ds.index.intersection(vaids)
        if len(tids) < MIN_TRAIN or len(vids) < MIN_VAL:
            log(
                "SKIP",
                f"{did} {ds.DRUG_NAME.iloc[0] if len(ds) else '?'}: TRAIN={len(tids)} VAL={len(vids)} below {MIN_TRAIN}/{MIN_VAL}",
            )
            continue
        ytr = ds.loc[tids, "AUC"].to_numpy(dtype="float32")
        yva = ds.loc[vids, "AUC"].to_numpy(dtype="float32")
        baseline = float(
            mean_absolute_error(yva, np.repeat(float(ytr.mean()), len(vids)))
        )
        tcan = cancer.reindex(tids)
        vcan = cancer.reindex(vids)
        cancer_means = (
            pd.DataFrame({"cancer": tcan.to_numpy(), "AUC": ytr})
            .groupby("cancer", dropna=True)
            .AUC.mean()
        )
        cancer_pred = vcan.map(cancer_means).fillna(float(ytr.mean())).to_numpy(float)
        cancer_mae = float(mean_absolute_error(yva, cancer_pred))
        for label, matrix in features.items():
            try:
                xtr = matrix.loc[tids].to_numpy(dtype="float32")
                xva = matrix.loc[vids].to_numpy(dtype="float32")
                mae, r2, pred, _ = trial_model(xtr, ytr, xva, yva, label)
                records.append(
                    dict(
                        drug_id=did,
                        drug_name=str(ds.DRUG_NAME.iloc[0]),
                        modality=label,
                        train_n=len(tids),
                        val_n=len(vids),
                        n_features=matrix.shape[1],
                        MAE=mae,
                        R2=r2,
                        MAE_train_mean=baseline,
                        MAE_cancer_mean=cancer_mae,
                        delta_vs_mean=baseline - mae,
                    )
                )
                predictions.extend(
                    dict(
                        model_id=id_,
                        drug_id=did,
                        drug_name=str(ds.DRUG_NAME.iloc[0]),
                        modality=label,
                        observed_auc=float(obs),
                        predicted_auc=float(p),
                    )
                    for id_, obs, p in zip(vids, yva, pred)
                )
            except Exception as e:
                log("ERROR", f"{did} {label}: {type(e).__name__}: {e}")
        log(
            "DRUG",
            f"{ds.DRUG_NAME.iloc[0]} [{did}]: matched TRAIN={len(tids)}, VAL={len(vids)}",
        )
    result = pd.DataFrame(records)
    if result.empty:
        raise RuntimeError("No models succeeded; inspect errors")
    result.to_csv(OUT / "multiomics_validation_comparison.csv", index=False)
    pd.DataFrame(predictions).to_csv(
        OUT / "multiomics_validation_predictions.csv", index=False
    )
    common = result.groupby("drug_id").modality.nunique()
    complete = common[common == len(features)].index
    paired = result[result.drug_id.isin(complete)]
    log(
        "RESULT",
        f"Successful records={len(result)}, drugs with all five approaches={len(complete)}; model comparisons are paired within same drug and sample cohort",
    )
    summary = (
        paired.groupby("modality")
        .agg(
            drugs=("drug_id", "nunique"),
            median_MAE=("MAE", "median"),
            median_R2=("R2", "median"),
            median_gain_vs_trainmean=("delta_vs_mean", "median"),
        )
        .reset_index()
        .sort_values("median_MAE")
    )
    summary.to_csv(OUT / "modality_summary.csv", index=False)
    log(
        "MODALITY SUMMARY",
        "\n"
        + fmt(
            summary,
            [
                "modality",
                "drugs",
                "median_MAE",
                "median_R2",
                "median_gain_vs_trainmean",
            ],
            10,
        ),
    )
    wide = paired.pivot(
        index=["drug_id", "drug_name"], columns="modality", values="MAE"
    ).reset_index()
    if {"AllFour", "RNA"}.issubset(wide.columns):
        wide["gain_all_vs_RNA"] = wide["RNA"] - wide["AllFour"]
        log(
            "COMPARISON",
            f"AllFour better than RNA in {int((wide.gain_all_vs_RNA > 0).sum())}/{len(wide)} drugs; median paired MAE improvement={wide.gain_all_vs_RNA.median():+.5f}",
        )
        log(
            "PER DRUG",
            "\n"
            + fmt(
                wide.sort_values("gain_all_vs_RNA", ascending=False),
                [
                    "drug_name",
                    "RNA",
                    "Protein",
                    "RNA+Protein",
                    "RNA+Mutation+CNV",
                    "AllFour",
                    "gain_all_vs_RNA",
                ],
                20,
            ),
        )
    wide.to_csv(OUT / "paired_drug_modality_MAE.csv", index=False)
    log(
        "BASELINES",
        f"Per-drug TRAIN-mean and cancer-type TRAIN-mean baselines included in multiomics_validation_comparison.csv; cancer types with no TRAIN representation use overall TRAIN mean",
    )
    log(
        "LIMITATIONS",
        "Validation subset has been reused for earlier exploratory work; no independent confirmatory inference. Protein intensities and CNV assays may carry source effects; missingness and coverage also limit inference.",
    )
    log(
        "LIMITATIONS",
        "Mutation summary predictors distinguish recorded events from no recorded event; NOT altered vs confirmed wild type. Unimputed protein and CNV missingness is quantified in saved feature matrix, median imputation fitted TRAIN only.",
    )
    log(
        "NEXT",
        "Inspect paired modality differences, per-drug gains, missingness and RNA/protein concordance; then targeted SHAP and confirmatory external PRISM evaluation where compound identity and endpoint are compatible.",
    )
    log("STATUS", f"SUCCESS | elapsed {(time.time() - T0) / 60:.1f} min | {OUT}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log("STATUS", f"FAILED: {type(e).__name__}: {e}")
        log("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 92 + "\n=== BEGIN PROJECT 5 PHASE 3A RESULTS ===")
        print("\n".join(REPORT))
        print("=== END PROJECT 5 PHASE 3A RESULTS ===\n" + "=" * 92)

# ----------------------------------------------------------------------------
# Explainability and modality contribution
# ----------------------------------------------------------------------------

from pathlib import Path
import traceback, time, warnings
import numpy as np
import pandas as pd
from xgboost import XGBRegressor, DMatrix
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import spearmanr

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
P2 = ROOT / "Results" / "Phase2"
P3A = ROOT / "Results" / "Phase3A"
OUT = ROOT / "Results" / "Phase3B"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261028
DRUG_IDS = {"1047": "Nutlin-3a (-)", "1003": "Camptothecin", "1060": "PD0325901"}
N_PERMUTATIONS = 12
N_STABILITY = 4
REPORT = []
START = time.time()


def log(tag, s):
    REPORT.append(f"[{tag}] {s}")


def norm(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def tab(df, n=15):
    return (
        df.head(n).to_string(index=False, float_format=lambda x: f"{x:.5g}")
        if len(df)
        else "None"
    )


# Phase 3A and 3B are consolidated in this file. Reuse the Phase 3A readers
# already defined above instead of importing a removed pre-consolidation file.
p3a = sys.modules[__name__]


def fit_model(xtrain, ytrain, seed):
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
    reg.fit(xtrain, ytrain)
    return reg


def preprocess(train, val):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Skipping features without any observed values"
        )
        imputer = SimpleImputer(strategy="median", keep_empty_features=True)
        a = imputer.fit_transform(train)
        b = imputer.transform(val)
    scaler = StandardScaler()
    return scaler.fit_transform(a).astype("float32"), scaler.transform(b).astype(
        "float32"
    )


def make_data():
    splits = pd.read_csv(P2 / "model_splits.csv")
    splits.model_id = norm(splits.model_id)
    if splits.model_id.duplicated().any():
        raise ValueError("Duplicate model_id in splits")
    rna = pd.read_csv(P2 / "expression_selected_genes.csv.gz", compression="gzip")
    rna.model_id = norm(rna.model_id)
    rna = rna.set_index("model_id").apply(pd.to_numeric, errors="coerce")
    ids = set(rna.index) & set(splits.model_id)
    pro = p3a.read_proteome(ids)
    cnv = p3a.read_cnv(ids)
    mut = p3a.read_mutation(ids)
    av = pd.read_csv(p3a.file("model_dataset_availability_"))
    av.model_id = norm(av.model_id)
    av = av.drop_duplicates("model_id").set_index("model_id")
    covcols = [
        c
        for c in av
        if c in ("Mutation Sanger Cell Lines (WES)", "Mutation Broad Cell Lines (WES)")
    ]
    if not covcols:
        raise ValueError("Mutation assay coverage flags not found")
    covered = pd.Series(False, index=av.index)
    for c in covcols:
        covered |= av[c].astype("string").str.lower().eq("true").fillna(False)
    matched = sorted(
        ids & set(pro.index) & set(cnv.index) & set(covered[covered].index)
    )
    s = splits.set_index("model_id").loc[matched, "split"]
    trids = s[s == "train"].index.tolist()
    vaids = s[s == "val"].index.tolist()
    if len(matched) != 814 or len(trids) != 569 or len(vaids) != 127:
        raise ValueError(
            f"Matched cohort differs from recorded Phase3A: matched={len(matched)} train={len(trids)} val={len(vaids)}. Reconcile before comparisons."
        )
    # Reuse the recorded Phase3A feature set verbatim; no re-selection using outcomes.
    feat = pd.read_csv(P3A / "selected_multiomics_features.csv")
    pc = feat.loc[feat.modality == "protein", "feature"].astype(str).tolist()
    cc = feat.loc[feat.modality == "CNV", "feature"].astype(str).tolist()
    mc = feat.loc[feat.modality == "mutation", "feature"].astype(str).tolist()
    if not pc or not cc or not mc:
        raise ValueError("Selected features missing from Phase3A")
    for name, mat, cols in [
        ("protein", pro, pc),
        ("CNV", cnv, cc),
        ("mutation", mut, mc),
    ]:
        missing = set(cols) - set(mat.columns)
        if missing:
            raise ValueError(
                f"{name}: {len(missing)} recorded Phase3A features missing from source"
            )
    mut = mut.reindex(index=matched, fill_value=0)
    rna = rna.reindex(matched).astype("float32")
    pro = pro.reindex(matched)[pc].astype("float32")
    cnv = cnv.reindex(matched)[cc].astype("float32")
    mut = mut[mc].astype("float32")
    allfour = pd.concat(
        [
            rna.add_prefix("RNA:"),
            pro.add_prefix("PRO:"),
            mut.add_prefix("MUT:"),
            cnv.add_prefix("CNV:"),
        ],
        axis=1,
    )
    rmc = pd.concat(
        [rna.add_prefix("RNA:"), mut.add_prefix("MUT:"), cnv.add_prefix("CNV:")], axis=1
    )
    features = {"RNA": rna, "RNA+Mutation+CNV": rmc, "AllFour": allfour}
    responses = pd.read_excel(
        p3a.file("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    responses["model_id"] = norm(responses.SANGER_MODEL_ID)
    responses["drug_id"] = norm(responses.DRUG_ID)
    responses["AUC"] = pd.to_numeric(responses.AUC, errors="coerce")
    responses = responses[
        responses.model_id.isin(matched) & responses.AUC.between(0, 1)
    ].copy()
    if responses.duplicated(["model_id", "drug_id"]).any():
        raise ValueError("Duplicate response per model/drug")
    return responses, s, trids, vaids, features


def shap_values(reg, x):
    # XGBoost TreeSHAP built in; no optional shap package required.
    contrib = reg.get_booster().predict(DMatrix(x), pred_contribs=True)
    margins = reg.get_booster().predict(DMatrix(x), output_margin=True)
    err = float(np.max(np.abs(contrib.sum(axis=1) - margins)))
    return contrib[:, :-1], contrib[:, -1], err


def main():
    log(
        "DESIGN",
        "Phase3A-preserved 814 matched models, fixed feature lists and XGBoost hyperparameters. TRAIN fits, VAL interpretation only; no TEST use.",
    )
    log(
        "DESIGN",
        "Three case studies selected from previously examined VAL performance; every inference is exploratory.",
    )
    log(
        "CAUTION",
        "SHAP values quantify contributions to model output, not causal or clinically validated biomarkers.",
    )
    responses, s, trids, vaids, features = make_data()
    log(
        "MATCH",
        f"models={len(s)} train={len(trids)} val={len(vaids)} test={(s == 'test').sum()}; modality dimensions={ {k: v.shape[1] for k, v in features.items()} }",
    )
    log(
        "FEATURE SOURCE",
        "Reused Results/Phase3A/selected_multiomics_features.csv; all RNA features from Phase2 original TRAIN-selected matrix",
    )
    metrics = []
    importance = []
    stability = []
    ablation = []
    assocs = []
    missingness = []
    predictions = []
    for did, name in DRUG_IDS.items():
        sub = responses[responses.drug_id == did].set_index("model_id")
        ti = [x for x in trids if x in sub.index]
        vi = [x for x in vaids if x in sub.index]
        if len(ti) < 100 or len(vi) < 20:
            raise ValueError(
                f"{name}: insufficient shared TRAIN / VAL: {len(ti)}, {len(vi)}"
            )
        ytr = sub.loc[ti, "AUC"].to_numpy("float32")
        yv = sub.loc[vi, "AUC"].to_numpy("float32")
        log(
            "DRUG",
            f"{name} [{did}] TRAIN={len(ti)}, VAL={len(vi)} AUC std VAL={np.std(yv):.5f}",
        )
        for label, frame in features.items():
            xtr_raw = frame.loc[ti]
            xv_raw = frame.loc[vi]
            for block in ["RNA", "PRO", "MUT", "CNV"]:
                cols = [
                    c
                    for c in frame
                    if (c.startswith(block + ":") if label != "RNA" else block == "RNA")
                ]
                if cols:
                    missingness.append(
                        dict(
                            drug_id=did,
                            drug_name=name,
                            model=label,
                            block=block,
                            feature_count=len(cols),
                            train_missing_rate=float(
                                xtr_raw[cols].isna().to_numpy().mean()
                            ),
                            val_missing_rate=float(
                                xv_raw[cols].isna().to_numpy().mean()
                            ),
                        )
                    )
            xtr, xv = preprocess(xtr_raw, xv_raw)
            reg = fit_model(xtr, ytr, SEED)
            yp = reg.predict(xv)
            mae = float(mean_absolute_error(yv, yp))
            r2 = float(r2_score(yv, yp))
            metrics.append(
                dict(
                    drug_id=did,
                    drug_name=name,
                    model=label,
                    train_n=len(ti),
                    val_n=len(vi),
                    MAE=mae,
                    R2=r2,
                )
            )
            predictions.extend(
                dict(
                    drug_id=did,
                    drug_name=name,
                    model=label,
                    model_id=i,
                    observed_auc=float(a),
                    predicted_auc=float(b),
                )
                for i, a, b in zip(vi, yv, yp)
            )
            if label not in ("AllFour", "RNA+Mutation+CNV"):
                continue
            shap, base, err = shap_values(reg, xv)
            if err > 1e-3:
                raise ValueError(f"{name}/{label} SHAP additivity deviation {err}")
            imp = np.mean(np.abs(shap), axis=0)
            for j, col in enumerate(frame.columns):
                importance.append(
                    dict(
                        drug_id=did,
                        drug_name=name,
                        model=label,
                        feature=col,
                        block=col.split(":", 1)[0],
                        mean_abs_SHAP=float(imp[j]),
                        mean_signed_SHAP=float(np.mean(shap[:, j])),
                    )
                )
            log("SHAP", f"{name}/{label}: max additivity margin error {err:.3g}")
            if label != "AllFour":
                continue
            # Correlation of raw molecular features with measured response: descriptive
            # only, particularly for missingness and mutation/no-record indicators.
            leaders = np.argsort(imp)[::-1][:25]
            for j in leaders:
                raw = xv_raw.iloc[:, j].to_numpy("float64")
                mask = np.isfinite(raw) & np.isfinite(yv)
                rho = (
                    float(spearmanr(raw[mask], yv[mask]).statistic)
                    if mask.sum() >= 12 and np.unique(raw[mask]).size > 2
                    else np.nan
                )
                assocs.append(
                    dict(
                        drug_id=did,
                        drug_name=name,
                        feature=frame.columns[j],
                        n=int(mask.sum()),
                        rho_auc=rho,
                        mean_abs_SHAP=float(imp[j]),
                    )
                )
            # Modality-block permutation ablation on VAL. Broken cross-feature
            # correlations make this model-reliance diagnostic, not a causal test.
            for block in ["RNA", "PRO", "MUT", "CNV"]:
                js = np.array(
                    [
                        j
                        for j, c in enumerate(frame.columns)
                        if c.startswith(block + ":")
                    ]
                )
                if not len(js):
                    continue
                diffs = []
                for p in range(N_PERMUTATIONS):
                    rng = np.random.default_rng(
                        SEED + int(did) * 31 + p * 7 + len(block)
                    )
                    order = rng.permutation(len(vi))
                    xx = xv.copy()
                    xx[:, js] = xv[order][:, js]
                    diffs.append(float(mean_absolute_error(yv, reg.predict(xx)) - mae))
                ablation.append(
                    dict(
                        drug_id=did,
                        drug_name=name,
                        block=block,
                        n_features=len(js),
                        delta_MAE_mean=float(np.mean(diffs)),
                        delta_MAE_sd=float(np.std(diffs, ddof=1)),
                        permutations=N_PERMUTATIONS,
                    )
                )
            # TRAIN resampling for feature-attribution stability. Do not run feature
            # selection on VAL, and do not use these samples for test-set inference.
            ranks = []
            for b in range(N_STABILITY):
                rng = np.random.default_rng(SEED + int(did) * 97 + b)
                sel = np.sort(
                    rng.choice(len(ti), size=int(0.8 * len(ti)), replace=False)
                )
                trsub = xtr_raw.iloc[sel]
                ysub = ytr[sel]
                a, v = preprocess(trsub, xv_raw)
                m = fit_model(a, ysub, SEED + b + 101)
                contribution, _, e = shap_values(m, v)
                if e > 1e-3:
                    raise ValueError("SHAP stability additivity failed")
                avg = np.mean(np.abs(contribution), axis=0)
                ranks.append(set(np.argsort(avg)[-20:]))
                for j in np.argsort(avg)[-20:]:
                    stability.append(
                        dict(
                            drug_id=did,
                            drug_name=name,
                            resample=b + 1,
                            feature=frame.columns[j],
                            mean_abs_SHAP=float(avg[j]),
                        )
                    )
            jacc = [
                len(ranks[i] & ranks[j]) / len(ranks[i] | ranks[j])
                for i in range(len(ranks))
                for j in range(i + 1, len(ranks))
            ]
            log(
                "STABILITY",
                f"{name}: 4 resamples, mean top20-feature Jaccard={np.mean(jacc):.3f} (higher=more reproducible)",
            )
        log(
            "DONE",
            f"{name}: completed paired baseline, SHAP, modality permutations, TRAIN subsample stability",
        )
    for fn, rows in [
        ("model_comparison.csv", metrics),
        ("shap_feature_importance.csv", importance),
        ("stability_top20_features.csv", stability),
        ("permutation_modality_ablation.csv", ablation),
        ("shap_top_feature_auc_associations.csv", assocs),
        ("feature_missingness.csv", missingness),
        ("validation_predictions.csv", predictions),
    ]:
        pd.DataFrame(rows).to_csv(OUT / fn, index=False)
    m = pd.DataFrame(metrics)
    imp = pd.DataFrame(importance)
    abl = pd.DataFrame(ablation)
    stab = pd.DataFrame(stability)
    log(
        "MODEL COMPARISON",
        "\n" + tab(m[["drug_name", "model", "MAE", "R2", "val_n"]], 20),
    )
    log(
        "TOP FEATURES",
        "\n"
        + tab(
            imp.sort_values(
                ["drug_name", "model", "mean_abs_SHAP"], ascending=[True, True, False]
            )
            .groupby(["drug_name", "model"], sort=False)
            .head(12)[["drug_name", "model", "feature", "mean_abs_SHAP"]],
            72,
        ),
    )
    log(
        "MODALITY ABLATION",
        "Positive delta MAE means permutation worsens prediction; negative values mean permutation unexpectedly improves it.\n"
        + tab(abl[["drug_name", "block", "delta_MAE_mean", "delta_MAE_sd"]], 30),
    )
    stable_counts = (
        stab.groupby(["drug_name", "feature"])["resample"]
        .nunique()
        .reset_index(name="resample_count")
    )
    log(
        "STABLE FEATURES",
        "\n"
        + tab(
            stable_counts[stable_counts.resample_count >= 3].sort_values(
                ["drug_name", "resample_count"], ascending=[True, False]
            ),
            45,
        ),
    )
    log(
        "LIMITATIONS",
        "SHAP values depend on the fitted model and correlated input features; not proof of molecular mechanism. Permutation ablation disturbs feature correlations.",
    )
    log(
        "LIMITATIONS",
        "Validation data reused in prior exploratory analyses; model and drug choices are not independently confirmed. RNA/protein/CNV source effects and mutation summary no-record ambiguity remain.",
    )
    log(
        "NEXT",
        "Use stable feature groups to formulate prespecified biological hypotheses, then independently evaluate compatible, chemically verified GDSC-PRISM drug pairs.",
    )
    log(
        "STATUS",
        f"SUCCESS; elapsed={(time.time() - START) / 60:.1f} min; outputs={OUT}",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log("STATUS", f"FAILED {type(e).__name__}: {e}")
        log("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 93 + "\n=== BEGIN PROJECT 5 PHASE 3B RESULTS ===")
        print("\n".join(REPORT))
        print("=== END PROJECT 5 PHASE 3B RESULTS ===\n" + "=" * 93)
