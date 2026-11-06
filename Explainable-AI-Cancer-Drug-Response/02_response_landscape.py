"""Quantify drug-response heterogeneity within and between cancer lineages.

Uses within-drug comparisons so fitted AUC values from different concentration ranges are not ranked as directly comparable potency measurements.

Run from the repository root with ``python 02_response_landscape.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# Cancer-lineage response heterogeneity
# ----------------------------------------------------------------------------

from pathlib import Path
import sys, traceback, warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase2C"
OUT.mkdir(parents=True, exist_ok=True)
MIN_CANCER_N = 8  # for individual cancer distribution summaries
MIN_TEST_GROUP_N = 10  # for between-cancer omnibus testing
MIN_GROUPS = 3
MIN_WITHIN_N = 12  # for relative sensitivity profiles
MIN_TEST_PRED_N = 8  # per-cancer held-out error reporting
TOP_N = 15
LOG = []


def say(tag, msg):
    LOG.append(f"[{tag}] {msg}")


def name_file(prefix):
    hits = list(DATA.glob(prefix + "*"))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one {prefix}*, got {[x.name for x in hits]}")
    return hits[0]


def norm(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def short_rows(df, cols, n=TOP_N, decimals=4):
    if df.empty:
        return "No eligible results."
    z = df.loc[:, cols].head(n).copy()
    for col in z.select_dtypes("number").columns:
        z[col] = z[col].round(decimals)
    return z.to_string(index=False)


def main():
    say(
        "DESIGN",
        "Observational cell-line drug-response screen; all compounds, no retrospective AI training.",
    )
    say(
        "DESIGN",
        f"Between-cancer test: >= {MIN_GROUPS} cancer groups with >= {MIN_TEST_GROUP_N} cell lines/group; BH-FDR across tested drugs.",
    )
    say(
        "DESIGN",
        "Within-cancer sensitive/less-sensitive labels are quartiles of observed AUC, relative only; NOT clinical sensitive/resistant calls.",
    )
    say(
        "DESIGN",
        "Lower GDSC AUC generally indicates stronger reduction in viability; compare cancers WITHIN a drug, not absolute AUC across drugs.",
    )
    model = pd.read_csv(
        name_file("model_list_"),
        usecols=[
            "model_id",
            "model_name",
            "patient_id",
            "tissue",
            "cancer_type",
            "model_type",
            "tissue_status",
        ],
        low_memory=False,
    )
    model["model_id"] = norm(model.model_id)
    if model.model_id.duplicated().any():
        raise ValueError("Duplicate model IDs in metadata")
    raw = pd.read_excel(name_file("GDSC2_fitted_dose_response_"))
    need = [
        "SANGER_MODEL_ID",
        "DRUG_ID",
        "DRUG_NAME",
        "AUC",
        "LN_IC50",
        "MIN_CONC",
        "MAX_CONC",
    ]
    missing = set(need) - set(raw.columns)
    if missing:
        raise ValueError(f"GDSC missing columns: {missing}")
    df = raw[need].copy()
    df["model_id"] = norm(df.SANGER_MODEL_ID)
    df["drug_id"] = norm(df.DRUG_ID)
    df["AUC"] = pd.to_numeric(df.AUC, errors="coerce")
    df["LN_IC50"] = pd.to_numeric(df.LN_IC50, errors="coerce")
    df = df.dropna(subset=["model_id", "drug_id", "AUC"])
    df = df[(df.AUC >= 0) & (df.AUC <= 1)].copy()
    df = df.merge(model, on="model_id", how="inner", validate="many_to_one")
    if df.duplicated(["drug_id", "model_id"]).any():
        raise ValueError(
            "Repeated drug-model observations found: review assay replicates instead of silently aggregating."
        )
    df["cancer_type"] = df.cancer_type.fillna("Unknown")
    df["DRUG_NAME"] = df.DRUG_NAME.fillna("Unnamed")
    say(
        "DATA",
        f"{len(df):,} fitted observations; {df.model_id.nunique():,} models; {df.drug_id.nunique():,} drugs; {df.cancer_type.nunique():,} cancer categories",
    )
    say(
        "DATA",
        f"Non-tumour models: {(df.cancer_type == 'Non-Cancerous').sum():,} response rows; excluded from tumour-only comparisons.",
    )
    # The primary cancer response analysis focuses on tumour-derived cell lines.
    df_t = df[
        (df.cancer_type != "Non-Cancerous")
        & (df.tissue_status.isin(["Tumour", "Metastasis"]))
        & (df.model_type == "Cell Line")
    ].copy()
    say(
        "COHORT",
        f"Tumour/metastasis cell-line observations: {len(df_t):,}; {df_t.model_id.nunique():,} models; {df_t.drug_id.nunique():,} drugs.",
    )
    df_t.to_csv(
        OUT / "gdsc2_tumour_cellline_drug_responses.csv.gz",
        index=False,
        compression="gzip",
    )

    per = (
        df_t.groupby(["drug_id", "DRUG_NAME", "cancer_type"], observed=True)
        .agg(
            n=("model_id", "nunique"),
            median_auc=("AUC", "median"),
            mean_auc=("AUC", "mean"),
            sd_auc=("AUC", "std"),
            q25=("AUC", lambda x: x.quantile(0.25)),
            q75=("AUC", lambda x: x.quantile(0.75)),
            min_auc=("AUC", "min"),
            max_auc=("AUC", "max"),
        )
        .reset_index()
    )
    per["iqr_auc"] = per.q75 - per.q25
    per.sort_values(["drug_id", "cancer_type"]).to_csv(
        OUT / "drug_by_cancer_auc_summary.csv", index=False
    )
    eligible = per[per.n >= MIN_CANCER_N]
    say(
        "CANCER",
        f"{len(per):,} compound × cancer combinations; {len(eligible):,} with n>={MIN_CANCER_N}; {eligible.drug_id.nunique():,} drugs represented.",
    )

    stats = []
    for did, d in df_t.groupby("drug_id", sort=False):
        ct = d.cancer_type.value_counts()
        groups = [
            d.loc[d.cancer_type.eq(k), "AUC"].to_numpy(dtype=float)
            for k in ct.index
            if ct[k] >= MIN_TEST_GROUP_N
        ]
        if len(groups) < MIN_GROUPS:
            continue
        # Degenerate all-tied groups must not receive a fabricated p value.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                H, p = kruskal(*groups)
        except ValueError:
            H, p = np.nan, np.nan
        N = sum(map(len, groups))
        k = len(groups)
        stats.append(
            {
                "drug_id": did,
                "drug_name": d.DRUG_NAME.iloc[0],
                "n_tested": N,
                "n_cancers": k,
                "H": H,
                "p": p,
                "epsilon2": max(0, (H - k + 1) / (N - k))
                if np.isfinite(H) and N > k
                else np.nan,
                "median_range_auc": float(
                    max(np.median(g) for g in groups)
                    - min(np.median(g) for g in groups)
                ),
            }
        )
    tests = pd.DataFrame(stats)
    if tests.empty:
        raise ValueError("No compounds eligible for between-cancer testing")
    tests["q_fdr"] = np.nan
    keep = tests.p.notna() & np.isfinite(tests.p)
    tests.loc[keep, "q_fdr"] = multipletests(tests.loc[keep, "p"], method="fdr_bh")[1]
    tests = tests.sort_values(["q_fdr", "epsilon2"], ascending=[True, False])
    tests.to_csv(OUT / "drug_between_cancer_kruskal_fdr.csv", index=False)
    say(
        "TESTS",
        f"{len(tests):,} drugs tested; {(tests.q_fdr < 0.05).sum():,} FDR<0.05; omnibus tests show differences, NOT which pairs differ.",
    )
    say(
        "TESTS",
        "Kruskal-Wallis significance may reflect different distribution shapes/spreads as well as medians; association, not causal effect.",
    )
    say(
        "TOP DIFFERENCES",
        "\n"
        + short_rows(
            tests,
            [
                "drug_name",
                "n_tested",
                "n_cancers",
                "median_range_auc",
                "epsilon2",
                "q_fdr",
            ],
        ),
    )
    top_effect = tests[(tests.q_fdr < 0.05) & (tests.n_cancers >= 4)].sort_values(
        "epsilon2", ascending=False
    )
    say(
        "TOP EFFECT SIZES",
        "\n"
        + short_rows(
            top_effect,
            [
                "drug_name",
                "n_tested",
                "n_cancers",
                "epsilon2",
                "median_range_auc",
                "q_fdr",
            ],
        ),
    )

    # Across-cancer patterns; retain only adequate group sample sizes; no pretending cell-line assays are clinical response.
    large = eligible[eligible.n >= MIN_TEST_GROUP_N]
    per_drug = (
        large.groupby(["drug_id", "DRUG_NAME"])
        .agg(
            n_cancers=("cancer_type", "size"),
            n_celllines=("n", "sum"),
            range_median=("median_auc", lambda v: v.max() - v.min()),
            min_median=("median_auc", "min"),
            max_median=("median_auc", "max"),
        )
        .reset_index()
    )
    per_drug = per_drug[per_drug.n_cancers >= MIN_GROUPS].sort_values(
        "range_median", ascending=False
    )
    per_drug.to_csv(OUT / "drug_across_cancer_median_ranges.csv", index=False)
    say(
        "TOP CANCER CONTRASTS",
        "\n"
        + short_rows(
            per_drug,
            [
                "DRUG_NAME",
                "n_cancers",
                "n_celllines",
                "range_median",
                "min_median",
                "max_median",
            ],
        ),
    )

    within = per[per.n >= MIN_WITHIN_N].sort_values("iqr_auc", ascending=False).copy()
    within.to_csv(OUT / "within_cancer_response_variability.csv", index=False)
    say(
        "WITHIN",
        f"{len(within):,} drug × cancer groups have n>={MIN_WITHIN_N}; IQR reflects observed heterogeneity, not assay reliability.",
    )
    say(
        "MOST VARIABLE WITHIN CANCER",
        "\n"
        + short_rows(
            within, ["DRUG_NAME", "cancer_type", "n", "median_auc", "iqr_auc"]
        ),
    )

    # Per drug × cancer: sample-specific relative quartile tags, only when n>=12.
    relative = []
    for (did, can), g in df_t.groupby(["drug_id", "cancer_type"], sort=False):
        if len(g) < MIN_WITHIN_N:
            continue
        q1, q3 = g.AUC.quantile([0.25, 0.75])
        if not q1 < q3:
            continue
        tagged = g[
            ["model_id", "model_name", "drug_id", "DRUG_NAME", "cancer_type", "AUC"]
        ].copy()
        tagged["relative_response"] = np.select(
            [tagged.AUC <= q1, tagged.AUC >= q3],
            ["lower_AUC_quartile", "higher_AUC_quartile"],
            default="middle_50_percent",
        )
        tagged["within_group_q25_auc"] = q1
        tagged["within_group_q75_auc"] = q3
        relative.append(tagged)
    rel = pd.concat(relative, ignore_index=True) if relative else pd.DataFrame()
    rel.to_csv(OUT / "within_cancer_relative_response_groups.csv", index=False)
    say(
        "CELL LINES",
        f"{len(rel):,} annotated cell-line measurements across {len(relative):,} eligible drug × cancer groups; no absolute sensitivity thresholds asserted.",
    )

    # Reuse existing held-out model predictions; no models are retrained and no test-driven feature selection.
    predict_file = ROOT / "Results" / "Phase2" / "heldout_test_predictions.csv"
    if predict_file.exists():
        pr = pd.read_csv(predict_file)
        needed = {"drug_name", "cancer_type", "observed_AUC", "predicted_AUC"}
        if not needed.issubset(pr):
            say(
                "PREDICTIONS",
                f"Skipped: heldout prediction file missing {needed - set(pr)}",
            )
        else:
            pr["abs_error"] = (pr.observed_AUC - pr.predicted_AUC).abs()
            pm = (
                pr.groupby(["drug_name", "cancer_type"])
                .agg(
                    n=("abs_error", "size"),
                    MAE=("abs_error", "mean"),
                    observed_median=("observed_AUC", "median"),
                    observed_IQR=(
                        "observed_AUC",
                        lambda s: s.quantile(0.75) - s.quantile(0.25),
                    ),
                )
                .reset_index()
            )
            pm = pm[pm.n >= MIN_TEST_PRED_N].sort_values(
                ["drug_name", "MAE"], ascending=[True, False]
            )
            pm.to_csv(OUT / "heldout_per_cancer_model_errors.csv", index=False)
            say(
                "PREDICTIONS",
                f"{len(pm):,} drug × cancer test groups with >= {MIN_TEST_PRED_N} held-out lines; per-group R2 not calculated for tiny groups.",
            )
            say(
                "HIGHEST HELDOUT ERRORS",
                "\n"
                + short_rows(
                    pm.sort_values("MAE", ascending=False),
                    ["drug_name", "cancer_type", "n", "MAE", "observed_IQR"],
                ),
            )
            say(
                "CAUTION",
                "Held-out test performance has already been inspected; this subgroup analysis is exploratory, not independent confirmation.",
            )
    else:
        say(
            "PREDICTIONS",
            "Previous Phase2 heldout_test_predictions.csv not found; biological analyses completed without it.",
        )

    # Specific examples of cross-cancer contrasts using observed medians (not drug-to-drug scores).
    examples = []
    for did in top_effect.head(10).drug_id:
        z = large[large.drug_id.eq(did)].sort_values("median_auc")
        if len(z) < 3:
            continue
        lo = z.iloc[0]
        hi = z.iloc[-1]
        examples.append(
            {
                "drug": lo.DRUG_NAME,
                "lower_auc_cancer": lo.cancer_type,
                "lower_auc": lo.median_auc,
                "lower_n": lo.n,
                "higher_auc_cancer": hi.cancer_type,
                "higher_auc": hi.median_auc,
                "higher_n": hi.n,
                "median_gap": hi.median_auc - lo.median_auc,
            }
        )
    ex = pd.DataFrame(examples)
    ex.to_csv(OUT / "example_cancer_median_contrasts.csv", index=False)
    say(
        "CONTRAST EXAMPLES",
        "\n"
        + (
            short_rows(ex, list(ex.columns), 10)
            if not ex.empty
            else "No contrasts eligible."
        ),
    )
    say(
        "LIMITS",
        "Screening-derived AUC differs from clinical drug efficacy; assays may vary by dose range and screen conditions.",
    )
    say(
        "LIMITS",
        "No biomarker hypothesis tests, post-hoc cancer pairs, independent PRISM comparisons or causal claims yet; these should follow prioritisation.",
    )
    say(
        "NEXT",
        "Use high-effect, sufficiently replicated drugs for within-cancer molecular association, and explicitly compare tissue/cancer effects with within-cancer effects.",
    )


if __name__ == "__main__":
    try:
        main()
        status = "SUCCESS"
    except Exception as e:
        status = "FAILED"
        say("ERROR", f"{type(e).__name__}: {e}")
        say("TRACEBACK", traceback.format_exc(limit=12))
        raise
    say("STATUS", status)
    (OUT / "phase2c_summary.txt").write_text("\n".join(LOG), encoding="utf-8")
    print("\n" + "=" * 88 + "\n=== BEGIN PROJECT 5 PHASE 2C RESULTS ===")
    print("\n".join(LOG))
    print("=== END PROJECT 5 PHASE 2C RESULTS ===\n" + "=" * 88)

# ----------------------------------------------------------------------------
# Within-lineage molecular associations
# ----------------------------------------------------------------------------

from pathlib import Path
import zipfile, io, traceback
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase2D"
OUT.mkdir(parents=True, exist_ok=True)
LOG = []
SEED = 20261028
# Drug name + cancer; explicitly resolve to unique DRUG_ID, do not combine same-name IDs.
TARGETS = [
    ("Dabrafenib", "Melanoma"),
    ("Docetaxel", "Breast Carcinoma"),
    ("Trametinib", "Melanoma"),
]
MIN_N = 20
MIN_FEATURE_N = 20
TOP = 12
# Hypothesis-driven candidate genes; explicit targeted readout, not claimed as validated biomarkers.
KNOWN = {
    "Dabrafenib": [
        "BRAF",
        "NRAS",
        "MAP2K1",
        "MAP2K2",
        "DUSP6",
        "SPRY2",
        "MITF",
        "PTEN",
    ],
    "Trametinib": [
        "BRAF",
        "NRAS",
        "KRAS",
        "MAP2K1",
        "MAP2K2",
        "DUSP6",
        "SPRY2",
        "PTEN",
    ],
    "Docetaxel": ["TUBB", "TUBB3", "ABCB1", "MAPT", "BCL2", "BCL2L1", "MKI67", "TOP2A"],
}


def say(section, msg):
    LOG.append(f"[{section}] {msg}")


def f(prefix):
    p = sorted(DATA.glob(prefix + "*"))
    if len(p) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {prefix}*: {[x.name for x in p]}"
        )
    return p[0]


def norm(v):
    return (
        v.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def qfmt(p):
    return "NA" if not np.isfinite(p) else f"{p:.3e}"


def show(df, cols, n=TOP):
    if df.empty:
        return "None"
    return (
        df.loc[:, cols]
        .head(n)
        .to_string(index=False, float_format=lambda v: f"{v:.4g}")
    )


def zip_read_selected(zip_path, match, **kw):
    with zipfile.ZipFile(zip_path) as z:
        cand = [x for x in z.namelist() if match in Path(x).name]
        if len(cand) != 1:
            raise RuntimeError(
                f"Cannot uniquely select {match} from {zip_path.name}: {cand}"
            )
        with z.open(cand[0]) as handle:
            return pd.read_csv(handle, **kw)


def get_expression(model_ids):
    zpath = f("rnaseq_merged_")
    member = "rnaseq_merged_rsem_tpm_"
    with zipfile.ZipFile(zpath) as z:
        members = [
            n
            for n in z.namelist()
            if member in Path(n).name and n.lower().endswith(".csv")
        ]
        if len(members) != 1:
            raise RuntimeError(f"Cannot select TPM matrix: {members}")
        with z.open(members[0]) as h:
            hdr = pd.read_csv(h, nrows=0).columns.tolist()
        actual = [x for x in hdr if x in model_ids]
        if not actual:
            raise RuntimeError("No matching expression model columns")
        select = hdr[:3] + actual
        with z.open(members[0]) as h:
            # Read only models in target analysis. The first three rows are metadata.
            expr = pd.read_csv(h, usecols=select, low_memory=False, dtype="string")
    if (
        expr.iloc[0, 0] != "model_name"
        or expr.iloc[1, 0] != "data_source"
        or expr.iloc[2, 0] != "gene_symbol"
    ):
        raise ValueError(
            "Expression metadata rows differ from Phase 1B; inspect schema"
        )
    source = expr.loc[1, actual].to_dict()
    genes = expr.iloc[3:].copy()
    genes = genes.rename(
        columns={
            hdr[0]: "gene_symbol",
            hdr[1]: "ensembl_gene_id",
            hdr[2]: "sanger_gene_id",
        }
    )
    genes["gene_symbol"] = genes["gene_symbol"].astype("string").str.strip()
    genes = genes[genes.gene_symbol.notna() & genes.gene_symbol.ne("")].copy()
    for c in actual:
        genes[c] = pd.to_numeric(genes[c], errors="coerce").astype("float32")
    # Most gene names unique. Drop duplicated names rather than silently pooling distinct genes.
    duplicated = int(genes.gene_symbol.duplicated(keep=False).sum())
    genes = genes.drop_duplicates("gene_symbol", keep=False).set_index("gene_symbol")[
        actual
    ]
    arr = genes.to_numpy(dtype=np.float32)
    arr[arr < 0] = np.nan
    genes.iloc[:, :] = np.log2(arr + 1)
    say(
        "RNA",
        f"TPM features={len(genes):,}; matching models={len(actual):,}; duplicated-symbol rows excluded={duplicated}; log2(TPM+1)",
    )
    return genes, source


def expr_associations(expr, df):
    cells = df.model_id.to_list()
    X = expr.loc[:, cells].to_numpy(dtype=np.float64)
    y = df.AUC.to_numpy(float)
    n = X.shape[1]
    valid_counts = np.isfinite(X).sum(axis=1)
    variances = np.nanvar(X, axis=1)
    usable = (valid_counts == n) & (variances > 1e-10)
    # To avoid modelling with distinct per-gene missingness, require all samples present.
    idx = np.flatnonzero(usable)
    x = X[idx]
    ranks = rankdata(x, axis=1, method="average")
    ry = rankdata(y, method="average")
    ranks -= ranks.mean(axis=1, keepdims=True)
    ry = ry - ry.mean()
    den = np.sqrt((ranks * ranks).sum(axis=1) * (ry * ry).sum())
    rho = np.divide(ranks @ ry, den, out=np.zeros(len(idx)), where=den > 0)
    # Standard asymptotic Spearman p-values; interpret small sample values cautiously.
    from scipy.stats import t as student_t

    tstat = rho * np.sqrt(np.maximum(n - 2, 1) / np.maximum(1 - rho * rho, 1e-12))
    p = np.minimum(1.0, 2 * student_t.sf(np.abs(tstat), df=n - 2))
    q = multipletests(p, method="fdr_bh")[1]
    ans = pd.DataFrame(
        {"gene": expr.index[idx].astype(str), "spearman_rho": rho, "p": p, "q_fdr": q}
    )
    ans["abs_rho"] = ans.spearman_rho.abs()
    ans = ans.sort_values(["q_fdr", "abs_rho"], ascending=[True, False])
    return ans


def mutations_for_models(model_ids):
    data = zip_read_selected(
        f("mutations_summary_"),
        "mutations_summary_",
        low_memory=False,
        usecols=["model_id", "gene_symbol", "cancer_driver", "effect"],
    )
    data["model_id"] = norm(data.model_id)
    return data[data.model_id.isin(model_ids)].copy()


def driver_associations(mut, df):
    # Presence is an annotated event; absence is not necessarily an independently verified negative.
    # Restrict mutation comparison to lines explicitly marked as covered by Sanger WES.
    av = pd.read_csv(
        f("model_dataset_availability_"),
        low_memory=False,
        usecols=["model_id", "Mutation Sanger Cell Lines (WES)"],
    )
    av.model_id = norm(av.model_id)
    ok = av.loc[
        av["Mutation Sanger Cell Lines (WES)"].astype(str).str.lower().eq("true"),
        "model_id",
    ]
    cov = df[df.model_id.isin(ok)].copy()
    if len(cov) < MIN_N:
        return pd.DataFrame(), len(cov)
    # NB: This is event-present vs no event found in the SUMMARY, not proven WT.
    observed = mut[mut.model_id.isin(cov.model_id)].dropna(subset=["gene_symbol"])
    counts = observed.groupby("gene_symbol").model_id.nunique()
    genes = counts[(counts >= 4) & (counts <= len(cov) - 4)].index
    out = []
    y = cov.set_index("model_id").AUC
    for gene in genes:
        pres = set(observed.loc[observed.gene_symbol.eq(gene), "model_id"])
        a = y[y.index.isin(pres)].to_numpy(float)
        b = y[~y.index.isin(pres)].to_numpy(float)
        if min(len(a), len(b)) < 4:
            continue
        pv = mannwhitneyu(a, b, alternative="two-sided").pvalue
        out.append(
            {
                "gene": gene,
                "event_n": len(a),
                "no_record_n": len(b),
                "median_auc_event": np.median(a),
                "median_auc_no_record": np.median(b),
                "delta_median_auc": np.median(a) - np.median(b),
                "p": pv,
            }
        )
    ans = pd.DataFrame(out)
    if not ans.empty:
        ans["q_fdr"] = multipletests(ans.p, method="fdr_bh")[1]
        ans = ans.sort_values(["q_fdr", "p"])
    return ans, len(cov)


def cnv_for_selected(model_ids, selected_genes):
    zpath = f("cnv_summary_")
    rows = []
    with zipfile.ZipFile(zpath) as z:
        n = [
            k
            for k in z.namelist()
            if "cnv_summary_" in Path(k).name and k.endswith(".csv")
        ]
        if len(n) != 1:
            raise RuntimeError("CNV member unavailable")
        with z.open(n[0]) as h:
            for chunk in pd.read_csv(
                h,
                usecols=["model_id", "symbol", "total_copy_number"],
                chunksize=180000,
                low_memory=False,
            ):
                chunk.model_id = norm(chunk.model_id)
                t = chunk[
                    chunk.model_id.isin(model_ids) & chunk.symbol.isin(selected_genes)
                ]
                if len(t):
                    rows.append(t)
    return (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(columns=["model_id", "symbol", "total_copy_number"])
    )


def cnv_associations(cnv, df):
    ans = []
    y = df.set_index("model_id").AUC
    for gene, g in cnv.groupby("symbol"):
        g = g[g.model_id.isin(y.index)].copy()
        g.total_copy_number = pd.to_numeric(g.total_copy_number, errors="coerce")
        # Do not silently select among mixed assays / repeated measurements.
        if g.model_id.duplicated().any():
            continue
        g = g.dropna(subset=["total_copy_number"])
        if len(g) < MIN_FEATURE_N or g.total_copy_number.nunique() < 4:
            continue
        yy = y.reindex(g.model_id).to_numpy(float)
        rho, p = spearmanr(g.total_copy_number.to_numpy(float), yy)
        if np.isfinite(p):
            ans.append({"gene": gene, "n": len(g), "rho": rho, "p": p})
    out = pd.DataFrame(ans)
    if not out.empty:
        out["q_fdr"] = multipletests(out.p, method="fdr_bh")[1]
        out = out.sort_values(["q_fdr", "p"])
    return out


def main():
    say(
        "DESIGN",
        "Exploratory within-cancer observational associations; NOT validation and NOT causal biomarker identification.",
    )
    say(
        "DESIGN",
        "Three preselected drug-cancer pairs from Phase 2C; no test-score-driven model refitting.",
    )
    say(
        "DESIGN",
        "Use DRUG_ID to disambiguate drugs; GDSC2 AUC comparisons are within individual drugs only.",
    )
    m = pd.read_csv(
        f("model_list_"),
        usecols=[
            "model_id",
            "model_name",
            "cancer_type",
            "tissue_status",
            "model_type",
        ],
        low_memory=False,
    )
    m.model_id = norm(m.model_id)
    if m.model_id.duplicated().any():
        raise ValueError("Duplicate model_id in annotation")
    drug = pd.read_excel(
        f("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    drug["model_id"] = norm(drug.SANGER_MODEL_ID)
    drug["drug_id"] = norm(drug.DRUG_ID)
    drug["AUC"] = pd.to_numeric(drug.AUC, errors="coerce")
    drug = drug.merge(m, on="model_id", how="inner", validate="many_to_one")
    drug = drug[
        (drug.model_type == "Cell Line")
        & drug.tissue_status.isin(["Tumour", "Metastasis"])
        & drug.AUC.between(0, 1)
    ].copy()
    if drug.duplicated(["model_id", "drug_id"]).any():
        raise ValueError("Repeated model-drug pairs")
    groups = []
    for d, c in TARGETS:
        tmp = drug[
            (drug.DRUG_NAME.str.casefold() == d.casefold()) & (drug.cancer_type == c)
        ]
        byid = tmp.groupby("drug_id").size().sort_values(ascending=False)
        if len(byid) != 1:
            say(
                "IDENTIFIER",
                f"{d} / {c}: found {len(byid)} drug IDs {byid.to_dict()}; choose no arbitrary version",
            )
            continue
        sub = tmp[tmp.drug_id == byid.index[0]].sort_values("model_id").copy()
        if len(sub) < MIN_N:
            say("COHORT", f"{d}/{c} skipped: n={len(sub)} <{MIN_N}")
            continue
        groups.append((d, c, byid.index[0], sub))
        say(
            "COHORT",
            f"{d} | {c} | DRUG_ID={byid.index[0]} | n={len(sub)} | median AUC={sub.AUC.median():.4f} | IQR={sub.AUC.quantile(0.75) - sub.AUC.quantile(0.25):.4f}",
        )
    if not groups:
        raise RuntimeError("No specified drug-cancer groups eligible")
    unique_models = set().union(*[set(x.model_id) for _, _, _, x in groups])
    say("COHORT", f"Unique cell lines across groups: {len(unique_models)}")
    expr, source = get_expression(unique_models)
    samples = []
    summary = []
    allgenes = []
    allmut = []
    allcnv = []
    muts = mutations_for_models(unique_models)
    cnv_genes = set.union(*[set(KNOWN[d]) for d, _, _, _ in groups])
    cnv = cnv_for_selected(unique_models, cnv_genes)
    for d, c, did, raw in groups:
        sub = raw[raw.model_id.isin(expr.columns)].copy().sort_values("model_id")
        say(
            "GROUP",
            f"{d} / {c}: expression-matched n={len(sub)}; Sanger/Broad RNA sources={pd.Series([source.get(i) for i in sub.model_id]).value_counts(dropna=False).to_dict()}",
        )
        if len(sub) < MIN_N:
            continue
        t1, t3 = sub.AUC.quantile([0.25, 0.75])
        sub["relative_response"] = np.where(
            sub.AUC <= t1,
            "lower_AUC_quartile",
            np.where(sub.AUC >= t3, "higher_AUC_quartile", "middle_50_percent"),
        )
        sub["drug_id"] = did
        sub["drug_name"] = d
        samples.append(
            sub[
                [
                    "drug_id",
                    "drug_name",
                    "cancer_type",
                    "model_id",
                    "model_name",
                    "AUC",
                    "relative_response",
                ]
            ]
        )
        es = expr_associations(expr, sub)
        es.insert(0, "drug_name", d)
        es.insert(1, "cancer_type", c)
        es.insert(2, "drug_id", did)
        allgenes.append(es)
        significant = int((es.q_fdr < 0.05).sum())
        say(
            "EXPRESSION",
            f"{d}/{c}: tested {len(es):,} variable complete genes; FDR<0.05: {significant}; Spearman exploratory, no adjustment for batch/source.",
        )
        say(
            "EXPRESSION TOP",
            f"{d}/{c}\n" + show(es, ["gene", "spearman_rho", "p", "q_fdr"]),
        )
        known = es[es.gene.isin(KNOWN[d])].copy()
        say(
            "TARGETED GENES",
            f"{d}/{c}\n" + show(known, ["gene", "spearman_rho", "p", "q_fdr"], n=20),
        )
        # This is gene-event presence, not a proven WT/mutant status comparison.
        ma, covered = driver_associations(muts, sub)
        if not ma.empty:
            ma.insert(0, "drug_name", d)
            ma.insert(1, "cancer_type", c)
            ma.insert(2, "drug_id", did)
            allmut.append(ma)
        say(
            "MUTATION",
            f"{d}/{c}: coverage-flagged WES lines={covered}; tested genes={len(ma)}; FDR<0.05={int((ma.q_fdr < 0.05).sum()) if len(ma) else 0}; event-present vs no-record ONLY.",
        )
        say(
            "MUTATION TOP",
            f"{d}/{c}\n"
            + (
                show(
                    ma,
                    [
                        "gene",
                        "event_n",
                        "no_record_n",
                        "delta_median_auc",
                        "p",
                        "q_fdr",
                    ],
                    n=10,
                )
                if len(ma)
                else "No eligible mutation contrasts"
            ),
        )
        ca = cnv_associations(cnv, sub)
        if not ca.empty:
            ca.insert(0, "drug_name", d)
            ca.insert(1, "cancer_type", c)
            ca.insert(2, "drug_id", did)
            allcnv.append(ca)
        say(
            "CNV",
            f"{d}/{c}: targeted genes tested={len(ca)}; FDR<0.05={int((ca.q_fdr < 0.05).sum()) if len(ca) else 0}; duplicates/mixed assay features excluded.",
        )
        say(
            "CNV TOP",
            f"{d}/{c}\n"
            + (
                show(ca, ["gene", "n", "rho", "p", "q_fdr"], n=12)
                if len(ca)
                else "No eligible CNV contrasts"
            ),
        )
        summary.append(
            {
                "drug": d,
                "drug_id": did,
                "cancer": c,
                "n": len(sub),
                "n_expression_genes": len(es),
                "n_expression_fdr05": significant,
                "n_mutation_genes": len(ma),
                "n_cnv_genes": len(ca),
            }
        )
    if samples:
        pd.concat(samples, ignore_index=True).to_csv(
            OUT / "selected_cellline_relative_response.csv", index=False
        )
    if allgenes:
        pd.concat(allgenes, ignore_index=True).to_csv(
            OUT / "within_cancer_expression_spearman_fdr.csv", index=False
        )
    if allmut:
        pd.concat(allmut, ignore_index=True).to_csv(
            OUT / "within_cancer_mutation_event_associations.csv", index=False
        )
    if allcnv:
        pd.concat(allcnv, ignore_index=True).to_csv(
            OUT / "within_cancer_selected_cnv_associations.csv", index=False
        )
    pd.DataFrame(summary).to_csv(OUT / "phase2d_group_summary.csv", index=False)
    say(
        "INTERPRETATION",
        "No clinical sensitivity classifications. Negative AUC association = higher expression associated with lower AUC, not causal sensitivity.",
    )
    say(
        "INTERPRETATION",
        "Mutation summary absent-event cannot establish verified wild type. Small n and transcriptomic source confounding require follow-up sensitivity analyses.",
    )
    say(
        "INTERPRETATION",
        "Genes/pathways nominated here are hypothesis generators; no formal pathway enrichment or independent PRISM validation performed.",
    )
    say(
        "NEXT",
        "Review statistically supported features and biological plausibility; select narrow hypotheses for matched multi-omics, XAI comparison and independent evaluation.",
    )


if __name__ == "__main__":
    try:
        main()
        status = "SUCCESS"
    except Exception as e:
        status = "FAILED"
        say("ERROR", f"{type(e).__name__}: {e}")
        say("TRACEBACK", traceback.format_exc(limit=12))
        raise
    say("STATUS", status)
    (OUT / "phase2d_summary.txt").write_text("\n".join(LOG), encoding="utf-8")
    print("\n" + "=" * 90 + "\n=== BEGIN PROJECT 5 PHASE 2D RESULTS ===")
    print("\n".join(LOG))
    print("=== END PROJECT 5 PHASE 2D RESULTS ===\n" + "=" * 90)
