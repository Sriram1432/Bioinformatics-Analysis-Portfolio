"""Investigate cancer-specific MAPK response patterns and predictive follow-up.

Combines relative response screens, a STRING-informed negative-control analysis and the selected six-gene feedback signature without claiming pathway activation.

Run from the repository root with ``python 07_mapk_followup.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# Focus-cancer comparative response screen
# ----------------------------------------------------------------------------

from pathlib import Path
import time, traceback, re
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase6A"
OUT.mkdir(parents=True, exist_ok=True)
P4 = ROOT / "Results" / "Phase4A"
P2 = ROOT / "Results" / "Phase2"
START = time.time()
REP = []
MIN_TARGET = 12
MIN_OTHER = 50
MIN_SHARED = 12
MIN_PRISM_CANCER = 8
# 'lung' is a reporting umbrella only: NSCLC/SCLC kept separate in tests.
GROUPS = {
    "Breast": "Breast Carcinoma",
    "Colorectal": "Colorectal Carcinoma",
    "Pancreatic": "Pancreatic Carcinoma",
    "NSCLC": "Non-Small Cell Lung Carcinoma",
    "SCLC": "Small Cell Lung Carcinoma",
}


def say(tag, msg):
    REP.append(f"[{tag}] {msg}")


def match_file(prefix):
    fs = list(DATA.glob(prefix + "*"))
    if len(fs) != 1:
        raise FileNotFoundError(
            f"Expected one Data/{prefix}*, found {[f.name for f in fs]}"
        )
    return fs[0]


def norm(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def col(df, names, need=True):
    cols = {re.sub("[^a-z0-9]", "", str(c).lower()): c for c in df.columns}
    for n in names:
        hit = cols.get(re.sub("[^a-z0-9]", "", n.lower()))
        if hit is not None:
            return hit
    if need:
        raise KeyError(f"Expected one of {names}; actual columns={list(df.columns)}")
    return None


def report(df, columns, n=15):
    if df.empty:
        return "No qualifying records"
    return df[columns].head(n).to_string(index=False, float_format=lambda x: f"{x:.5g}")


def correction(df, p_col="p", q_col="q_fdr"):
    df[q_col] = np.nan
    valid = df[p_col].notna()
    if valid.any():
        df.loc[valid, q_col] = multipletests(df.loc[valid, p_col], method="fdr_bh")[1]
    return df


def main():
    say(
        "DESIGN",
        "Tumour-derived CELL LINES only. NSCLC and SCLC analysed separately; no pooling of biologically different lung subtypes.",
    )
    say(
        "DESIGN",
        "All 295 GDSC2 drug IDs screened; within-drug cancer versus other-cancers rank effects. Lower AUC suggests stronger screened viability effect, not clinical efficacy.",
    )
    say(
        "DESIGN",
        "Shared candidates require evidence in >=2 of four focus organs (lung subtypes reported separately). BH correction across all evaluated focus contrasts.",
    )
    say(
        "DESIGN",
        "Drug AUCs are NOT compared numerically across compounds; sample-size and absolute within-drug distributions saved for inspection.",
    )
    metadata = pd.read_csv(match_file("model_list_"), low_memory=False)
    mid = col(metadata, ["model_id"])
    cancer = col(metadata, ["cancer_type"])
    mt = col(metadata, ["model_type"])
    ts = col(metadata, ["tissue_status"])
    mapping = col(
        metadata, ["broad_id", "depmap_id", "BROAD_ID", "DepMap_ID"], need=False
    )
    m = metadata[
        [mid, cancer, mt, ts] + ([mapping] if mapping is not None else [])
    ].copy()
    m = m.rename(
        columns={
            mid: "model_id",
            cancer: "cancer_type",
            mt: "model_type",
            ts: "tissue_status",
            **({mapping: "depmap_id"} if mapping is not None else {}),
        }
    )
    m["model_id"] = norm(m.model_id)
    if "depmap_id" in m:
        m["depmap_id"] = norm(m.depmap_id)
    m = m.dropna(subset=["model_id"]).drop_duplicates()
    if m.model_id.duplicated().any():
        raise ValueError("Non-unique model metadata; inspect model_list file")
    raw = pd.read_excel(
        match_file("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    raw["model_id"] = norm(raw.SANGER_MODEL_ID)
    raw["drug_id"] = norm(raw.DRUG_ID)
    raw["auc_gdsc"] = pd.to_numeric(raw.AUC, errors="coerce")
    raw = raw[raw.auc_gdsc.between(0, 1) & raw.model_id.notna()].merge(
        m, on="model_id", how="inner", validate="many_to_one"
    )
    raw = raw[
        (raw.model_type == "Cell Line")
        & raw.tissue_status.isin(["Tumour", "Metastasis"])
        & raw.cancer_type.notna()
    ].copy()
    if raw.duplicated(["model_id", "drug_id"]).any():
        raise ValueError(
            "Duplicate GDSC model-drug observations; do not treat as independent"
        )
    raw["focus"] = raw.cancer_type.map({v: k for k, v in GROUPS.items()})
    focus = raw[raw.focus.notna()].copy()
    counts = (
        m[
            (m.model_type == "Cell Line")
            & m.tissue_status.isin(["Tumour", "Metastasis"])
        ]
        .groupby("cancer_type")
        .model_id.nunique()
        .rename("models_metadata")
        .reset_index()
    )
    observed = (
        raw.groupby("cancer_type")
        .model_id.nunique()
        .rename("models_with_drug_response")
        .reset_index()
    )
    counts = counts.merge(observed, on="cancer_type", how="left")
    counts.to_csv(OUT / "all_cancer_cell_line_counts.csv", index=False)
    focus_counts = (
        focus.groupby(["focus", "cancer_type"])
        .agg(
            models=("model_id", "nunique"),
            measurements=("auc_gdsc", "size"),
            drugs=("drug_id", "nunique"),
        )
        .reset_index()
    )
    focus_counts.to_csv(OUT / "focus_cancer_coverage.csv", index=False)
    say(
        "COVERAGE",
        f"All tumour data: {len(raw):,} measurements, {raw.model_id.nunique()} cell lines, {raw.drug_id.nunique()} drug IDs, {raw.cancer_type.nunique()} cancer types.",
    )
    say(
        "FOCUS",
        "\n"
        + report(
            focus_counts,
            ["focus", "cancer_type", "models", "measurements", "drugs"],
            10,
        ),
    )
    # Cancer-vs-rest enrichment for every focus subtype and each drug.
    rows = []
    summary = []
    for did, drug in raw.groupby("drug_id", sort=False):
        nm = str(drug.DRUG_NAME.iloc[0])
        all_n = len(drug)
        for short, cat in GROUPS.items():
            sub = drug[drug.cancer_type == cat].auc_gdsc.to_numpy(float)
            other = drug[drug.cancer_type != cat].auc_gdsc.to_numpy(float)
            if len(sub) < MIN_TARGET or len(other) < MIN_OTHER:
                continue
            med = float(np.median(sub))
            med_other = float(np.median(other))
            u, p = mannwhitneyu(
                sub, other, alternative="two-sided", method="asymptotic"
            )
            prob = float(1 - u / (len(sub) * len(other)))
            effect = 2 * prob - 1
            rows.append(
                dict(
                    drug_id=did,
                    drug_name=nm,
                    focus=short,
                    cancer_type=cat,
                    n_cancer=len(sub),
                    n_other=len(other),
                    cancer_median_AUC=med,
                    cancer_AUC_IQR=float(np.subtract(*np.percentile(sub, [75, 25]))),
                    other_median_AUC=med_other,
                    median_gap=med - med_other,
                    probability_lower_AUC=prob,
                    rank_effect=effect,
                    p=float(p),
                )
            )
    enrich = correction(pd.DataFrame(rows))
    if enrich.empty:
        raise RuntimeError("No eligible focus cancer comparisons")
    enrich = enrich.sort_values(["q_fdr", "rank_effect"], ascending=[True, False])
    enrich.to_csv(OUT / "focus_drug_cancer_enrichment_all.csv", index=False)
    strong = enrich[
        (enrich.q_fdr < 0.05) & (enrich.rank_effect >= 0.3) & (enrich.median_gap < 0)
    ].copy()
    strong.to_csv(OUT / "focus_lower_auc_enrichment.csv", index=False)
    say(
        "TESTS",
        f"{len(enrich)} focus cancer-versus-other comparisons; significant lower AUC at FDR<.05 and rank_effect>=.30: {len(strong)}",
    )
    say(
        "TOP BY CANCER",
        "\n"
        + report(
            strong.sort_values("rank_effect", ascending=False),
            [
                "focus",
                "drug_name",
                "drug_id",
                "n_cancer",
                "cancer_median_AUC",
                "other_median_AUC",
                "rank_effect",
                "q_fdr",
            ],
            20,
        ),
    )
    for cat in GROUPS:
        t = strong[strong.focus.eq(cat)].sort_values(
            ["rank_effect", "n_cancer"], ascending=False
        )
        say(
            "FOCUS " + cat,
            "\n"
            + report(
                t,
                [
                    "drug_name",
                    "drug_id",
                    "n_cancer",
                    "cancer_median_AUC",
                    "other_median_AUC",
                    "rank_effect",
                    "q_fdr",
                ],
                12,
            ),
        )
    # Shared activity patterns: require per-cancer *independent comparisons against remaining cancer lines*, not sum of unlike drug AUCs.
    per = enrich.copy()
    per["qualified"] = (
        (per.q_fdr < 0.05) & (per.rank_effect >= 0.3) & (per.median_gap < 0)
    )
    per["organ"] = per.focus.map(
        {
            "NSCLC": "Lung",
            "SCLC": "Lung",
            "Breast": "Breast",
            "Colorectal": "Colorectal",
            "Pancreatic": "Pancreatic",
        }
    )
    # Lung counts once; choose strongest qualified subtype and preserve subtype in report.
    passing = (
        per[per.qualified]
        .sort_values("rank_effect", ascending=False)
        .drop_duplicates(["drug_id", "organ"])
    )
    summary = passing.groupby(["drug_id", "drug_name"], as_index=False).agg(
        n_organs=("organ", "nunique"),
        organs=("organ", lambda x: "; ".join(sorted(x))),
        subtype_details=("focus", lambda x: "; ".join(x)),
        min_rank_effect=("rank_effect", "min"),
        median_rank_effect=("rank_effect", "median"),
        min_organ_n=("n_cancer", "min"),
        max_q=("q_fdr", "max"),
    )
    summary = summary[summary.n_organs >= 2].sort_values(
        ["n_organs", "min_rank_effect", "min_organ_n"], ascending=False
    )
    summary.to_csv(OUT / "shared_lower_auc_compound_leads.csv", index=False)
    pan = summary[summary.organs.str.contains("Pancreatic", regex=False)]
    say(
        "SHARED",
        f"{len(summary)} drug IDs exhibit lower-AUC enrichment in >=2 focus organs; {len(pan)} involve pancreatic carcinoma.",
    )
    say(
        "SHARED TOP",
        "\n"
        + report(
            summary,
            [
                "drug_name",
                "drug_id",
                "n_organs",
                "organs",
                "min_rank_effect",
                "min_organ_n",
                "max_q",
            ],
            25,
        ),
    )
    say(
        "PANCREATIC SHARED",
        "\n"
        + report(
            pan,
            ["drug_name", "drug_id", "n_organs", "organs", "min_rank_effect", "max_q"],
            20,
        ),
    )
    # Cancer-by-drug standardized within-drug signature: NOT a between-drug ranking of effect potency.
    cr = (
        focus.groupby(["drug_id", "DRUG_NAME", "focus"])
        .auc_gdsc.agg(
            n="size",
            median_AUC="median",
            IQR_AUC=lambda x: x.quantile(0.75) - x.quantile(0.25),
        )
        .reset_index()
    )
    cr = cr.rename(columns={"DRUG_NAME": "drug_name"}).merge(
        enrich[["drug_id", "focus", "rank_effect", "q_fdr", "other_median_AUC"]],
        on=["drug_id", "focus"],
        how="left",
        validate="one_to_one",
    )
    cr.to_csv(OUT / "focus_cancer_drug_response_matrix_long.csv", index=False)
    cr.pivot(
        index=["drug_id", "drug_name"], columns="focus", values="median_AUC"
    ).to_csv(OUT / "focus_cancer_median_AUC_matrix.csv")
    # Inspect similarity between focus cancers only across identical drugs, with centred within-drug medians.
    med = cr[cr.n >= MIN_TARGET].pivot(
        index="drug_id", columns="focus", values="median_AUC"
    )
    centered = med.sub(med.median(axis=1), axis=0)
    simrows = []
    for i, a in enumerate(GROUPS):
        for b in list(GROUPS)[i + 1 :]:
            both = centered[[a, b]].dropna()
            rho = (
                float(spearmanr(both[a], both[b]).statistic)
                if len(both) >= 30 and both[a].nunique() > 1 and both[b].nunique() > 1
                else np.nan
            )
            simrows.append(
                dict(cancer_a=a, cancer_b=b, n_shared_drugs=len(both), rho=rho)
            )
    sim = pd.DataFrame(simrows).sort_values("rho", ascending=False)
    sim.to_csv(OUT / "focus_pairwise_drug_response_correlations.csv", index=False)
    say(
        "SIMILARITY",
        "\n" + report(sim, ["cancer_a", "cancer_b", "n_shared_drugs", "rho"], 12),
    )
    # Cancer-specific PRISM concordance on SAME models: only Phase4A provisional pairs, no chemical certification.
    cp = P4 / "candidate_pairs_for_analysis.csv"
    if cp.exists() and "depmap_id" in raw.columns:
        candidates = pd.read_csv(cp, dtype={"drug_id": str, "broad_id": str})
        candidates["drug_id"] = norm(candidates.drug_id)
        candidates["broad_id"] = norm(candidates.broad_id)
        ids = set(candidates.broad_id.dropna())
        source = match_file("secondary-screen-dose-response-curve-parameters")
        fields = pd.read_csv(source, nrows=0).columns.tolist()
        use = ["broad_id", "depmap_id", "auc"] + [
            v for v in ("passed_str_profiling", "screen_id") if v in fields
        ]
        pieces = []
        for chunk in pd.read_csv(
            source, usecols=use, chunksize=125000, low_memory=False
        ):
            chunk["broad_id"] = norm(chunk.broad_id)
            chunk = chunk[chunk.broad_id.isin(ids)].copy()
            if chunk.empty:
                continue
            chunk["depmap_id"] = norm(chunk.depmap_id)
            chunk["prism_auc"] = pd.to_numeric(chunk.auc, errors="coerce")
            chunk = chunk[chunk.prism_auc.notna() & chunk.depmap_id.notna()]
            if "passed_str_profiling" in chunk:
                ok = (
                    chunk.passed_str_profiling.astype("string")
                    .str.casefold()
                    .eq("true")
                    .fillna(False)
                )
                chunk = chunk[ok]
            pieces.append(chunk[["broad_id", "depmap_id", "prism_auc"]])
        if pieces:
            pp = pd.concat(pieces, ignore_index=True)
            duplicate = (pp.groupby(["broad_id", "depmap_id"]).size() > 1).sum()
            pp = pp.groupby(
                ["broad_id", "depmap_id"], as_index=False
            ).prism_auc.median()
            g = (
                focus[["model_id", "depmap_id", "drug_id", "focus", "auc_gdsc"]]
                .dropna(subset=["depmap_id"])
                .copy()
            )
            if g.duplicated(["drug_id", "depmap_id"]).any():
                ambiguous = g.groupby(["drug_id", "depmap_id"]).model_id.nunique()
                bad = ambiguous[ambiguous > 1].index
                g = g[~pd.MultiIndex.from_frame(g[["drug_id", "depmap_id"]]).isin(bad)]
                g = g.drop_duplicates(["drug_id", "depmap_id"])
            rows = []
            for pair in (
                candidates[["drug_id", "broad_id"]]
                .drop_duplicates()
                .itertuples(index=False)
            ):
                x = g[g.drug_id == pair.drug_id].merge(
                    pp[pp.broad_id == pair.broad_id],
                    on="depmap_id",
                    validate="one_to_one",
                )
                for name, z in x.groupby("focus"):
                    if len(z) < MIN_PRISM_CANCER:
                        continue
                    rho, p = (
                        spearmanr(z.auc_gdsc, z.prism_auc)
                        if z.auc_gdsc.nunique() > 1 and z.prism_auc.nunique() > 1
                        else (np.nan, np.nan)
                    )
                    rows.append(
                        dict(
                            drug_id=pair.drug_id,
                            broad_id=pair.broad_id,
                            focus=name,
                            n_shared=len(z),
                            rho=float(rho),
                            p=float(p),
                        )
                    )
            rc = (
                correction(pd.DataFrame(rows))
                if rows
                else pd.DataFrame(
                    columns=[
                        "drug_id",
                        "broad_id",
                        "focus",
                        "n_shared",
                        "rho",
                        "p",
                        "q_fdr",
                    ]
                )
            )
            rc.to_csv(
                OUT / "PRISM_within_cancer_cross_platform_concordance_PROVISIONAL.csv",
                index=False,
            )
            say(
                "PRISM",
                f"{len(rc)} matched drug-cancer assay comparisons with >= {MIN_PRISM_CANCER} lines; {duplicate} PRISM compound×cell-line combinations had >1 assay entries (median aggregation, exploratory).",
            )
            say(
                "PRISM TOP",
                "\n"
                + report(
                    rc[rc.rho > 0].sort_values("rho", ascending=False),
                    ["drug_id", "focus", "n_shared", "rho", "q_fdr"],
                    15,
                ),
            )
            # Add experimental concordance to shared leads without implying drug-cancer replication.
            if not rc.empty:
                joined = summary.merge(
                    rc[["drug_id", "focus", "rho", "n_shared", "broad_id"]],
                    on="drug_id",
                    how="left",
                )
                joined.to_csv(
                    OUT / "shared_leads_with_PRISM_subtype_evidence.csv", index=False
                )
    else:
        say(
            "PRISM",
            "Phase4A candidate map or DepMap cell IDs unavailable; PRISM follow-up skipped, not replaced with name guesses.",
        )
    # Prioritisation is based on differences, not evidence of therapeutic drug action.
    dd = raw[["drug_id", "DRUG_NAME"]].drop_duplicates()
    amb = dd.groupby("DRUG_NAME").drug_id.nunique()
    say(
        "IDENTIFIERS",
        f"{sum(amb > 1)} compound names correspond to multiple GDSC drug IDs: always keep drug_id.",
    )
    say(
        "PPI NEXT",
        "Test network-based features against RNA/proteomics baseline on fixed matched cell lines only AFTER obtaining a versioned high-confidence human PPI edge list; randomised-network control mandatory.",
    )
    say(
        "LIMITS",
        "Within-drug AUC enrichment ≠ absolute potency, efficacy, drug-selectivity or established repurposing. Screen concentrations differ across drugs. FDR-tested contrasts remain observational and selected.",
    )
    say(
        "LIMITS",
        "PRISM candidates remain name-matched only; identical cell lines across platforms, heterogeneous screen protocols, and replicate aggregation preclude independent generalisation claims.",
    )
    say("STATUS", f"SUCCESS | elapsed={time.time() - START:.1f}s | output={OUT}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        say("ERROR", f"{type(exc).__name__}: {exc}")
        say("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 108)
        print("=== BEGIN PROJECT 5 PHASE 6A RESULTS ===")
        for s in REP:
            print(s)
        print("=== END PROJECT 5 PHASE 6A RESULTS ===")
        print("=" * 108)

# ----------------------------------------------------------------------------
# Corrected PPI-informed negative-control analysis
# ----------------------------------------------------------------------------

from pathlib import Path
import importlib.util, urllib.request, gzip, time, traceback, warnings, json, re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy import sparse
from sklearn.model_selection import GroupKFold
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error
from xgboost import XGBRegressor

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
P2 = ROOT / "Results" / "Phase2"
OUT = ROOT / "Results" / "Phase6B"
PPI = DATA / "PPI"
OUT.mkdir(parents=True, exist_ok=True)
PPI.mkdir(parents=True, exist_ok=True)
SEED = 20261028
REPEATS = 2
FOLDS = 4
CUTOFF = 700
DRUGS = {"1060": "PD0325901", "1564": "SCH772984", "1014": "Refametinib"}
TARGETS = ["Colorectal Carcinoma", "Pancreatic Carcinoma"]
GROUPS = {
    "Colorectal Carcinoma": "Colorectal",
    "Pancreatic Carcinoma": "Pancreatic",
    "Breast Carcinoma": "Breast",
    "Non-Small Cell Lung Carcinoma": "NSCLC",
    "Small Cell Lung Carcinoma": "SCLC",
}
URL_ROOT = "https://stringdb-static.org/download"
FILE_INFO = "9606.protein.info.v12.5.txt.gz"
FILE_LINKS = "9606.protein.physical.links.v12.5.txt.gz"
LOG = []
START = time.time()


def say(k, v):
    LOG.append(f"[{k}] {v}")


def norm(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def fetch(filename, category):
    dst = PPI / filename
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    url = f"{URL_ROOT}/{category}/{filename}"
    try:
        urllib.request.urlretrieve(url, str(dst))
        with gzip.open(dst, "rb") as fh:
            fh.read(256)
    except Exception as e:
        dst.unlink(missing_ok=True)
        raise RuntimeError(
            f"Cannot download {url}: {e}. Download the file from https://string-db.org/cgi/download?species=Homo+sapiens and place it in Data/PPI/"
        )
    return dst


def helpers():
    fn = ROOT / "06_nested_validation.py"
    if not fn.exists():
        raise FileNotFoundError(f"Required earlier script missing: {fn}")
    spec = importlib.util.spec_from_file_location("phase5b_p5_helper", fn)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_ppi():
    info = fetch(FILE_INFO, "protein.info.v12.5")
    links = fetch(FILE_LINKS, "protein.physical.links.v12.5")
    inf = pd.read_csv(info, sep="\t", compression="gzip", low_memory=False)
    idcol = next((c for c in inf if "string_protein_id" in c.lower()), None)
    symcol = next((c for c in inf if "preferred_name" in c.lower()), None)
    if not idcol or not symcol:
        raise ValueError(f"Unexpected STRING protein.info columns: {list(inf.columns)}")
    inf["symbol"] = inf[symcol].astype(str).str.upper().str.strip()
    inf = inf[inf.symbol.str.match(r"^[A-Z0-9][A-Z0-9_.-]*$", na=False)]
    # Avoid protein-ID to symbol multi-mappings; collapse only unambiguous protein records.
    mapper = inf.drop_duplicates(idcol).set_index(idcol).symbol.to_dict()
    edges = []
    raw_n = 0
    for chunk in pd.read_csv(links, sep=r"\s+", compression="gzip", chunksize=250000):
        if not {"protein1", "protein2", "combined_score"}.issubset(chunk.columns):
            raise ValueError(
                f"Unexpected STRING physical links columns: {list(chunk.columns)}"
            )
        chunk = chunk[chunk.combined_score >= CUTOFF]
        raw_n += len(chunk)
        if chunk.empty:
            continue
        a = chunk.protein1.map(mapper)
        b = chunk.protein2.map(mapper)
        ok = a.notna() & b.notna() & a.ne(b)
        if ok.any():
            edges.append(
                pd.DataFrame(
                    {
                        "a": a[ok].to_numpy(),
                        "b": b[ok].to_numpy(),
                        "score": chunk.loc[ok, "combined_score"].to_numpy(
                            dtype=np.int16
                        ),
                    }
                )
            )
    if not edges:
        raise RuntimeError("No STRING links passing cutoff; inspect download")
    edge = pd.concat(edges, ignore_index=True)
    # Deduplicate undirected edges (possibly multiple proteins map to same symbol)
    a = edge.a.to_numpy(dtype=str)
    b = edge.b.to_numpy(dtype=str)
    x = np.where(a <= b, a, b)
    y = np.where(a <= b, b, a)
    edge = pd.DataFrame({"a": x, "b": y, "score": edge.score.to_numpy()})
    edge = edge.groupby(["a", "b"], as_index=False).score.max()
    edge.to_csv(
        OUT / "STRING_high_confidence_physical_edges_gene_symbol.csv.gz",
        index=False,
        compression="gzip",
    )
    say(
        "PPI",
        f"STRING v12.5 Homo sapiens physical links score>={CUTOFF}: {len(edge):,} unique gene-symbol edges, {len(set(edge.a) | set(edge.b)):,} genes. Evidence may include indirect physical associations.",
    )
    return edge


def network_matrix(columns, edge, seed=None):
    # Edges restricted to measured fold-selected predictors; no target-based edge selection.
    columns = list(columns)
    pos = {k: i for i, k in enumerate(columns)}
    e = edge[edge.a.isin(pos) & edge.b.isin(pos)]
    if seed is not None:
        # Label-permuted graph: same topology as retained graph, labels shuffled. No false claim of gene-specific degree preservation.
        rng = np.random.default_rng(seed)
        shuffled = np.asarray(columns, dtype=object).copy()
        rng.shuffle(shuffled)
        relabel = dict(zip(columns, shuffled))
        e = e.assign(a=e.a.map(relabel), b=e.b.map(relabel))
    i = e.a.map(pos).to_numpy(dtype=int)
    j = e.b.map(pos).to_numpy(dtype=int)
    val = e.score.to_numpy(dtype=float) / 1000.0
    n = len(columns)
    adj = sparse.csr_matrix((np.r_[val, val], (np.r_[i, j], np.r_[j, i])), shape=(n, n))
    degree = np.asarray(adj.sum(axis=1)).ravel()
    W = sparse.diags(1 / np.maximum(degree, 1e-9)) @ adj
    return W, int((degree > 0).sum()), len(e)


def model(seed):
    return XGBRegressor(
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


def main():
    say(
        "DESIGN",
        "Original TRAIN only; 2x4 patient-grouped CV. RNA/protein/CNV/mutation features selected inside folds; imputation fitted to training folds.",
    )
    say(
        "DESIGN",
        "Three selected MAPK drugs (not all 295); same cell lines for matched baseline, true PPI and shuffled-label control. Cancer effects evaluated on OOF predictions.",
    )
    say(
        "DESIGN",
        "PPI neighbourhood features summarize other measured molecules interacting in STRING physical network. Transductive external network only, no response-based edges.",
    )
    ppi = load_ppi()
    h = helpers()
    splits = pd.read_csv(P2 / "model_splits.csv", dtype={"model_id": "string"})
    splits.model_id = norm(splits.model_id)
    if splits.model_id.duplicated().any():
        raise ValueError("Nonunique split model IDs")
    if "group_key" not in splits:
        raise KeyError("Patient group key missing")
    original = set(splits.loc[splits.split.eq("train"), "model_id"])
    rna = h.read_all_rna(original)
    mod = h.phase3()
    pro = mod.read_proteome(original)
    cnv = mod.read_cnv(original)
    mut = mod.read_mutation(original)
    av = pd.read_csv(h.datafile("model_dataset_availability_"), low_memory=False)
    av.model_id = norm(av.model_id)
    av = av.drop_duplicates("model_id").set_index("model_id")
    covcols = [
        c
        for c in ("Mutation Sanger Cell Lines (WES)", "Mutation Broad Cell Lines (WES)")
        if c in av
    ]
    if not covcols:
        raise ValueError("Missing mutation WES coverage metadata")
    covered = pd.Series(False, index=av.index)
    for c in covcols:
        covered |= av[c].astype("string").str.casefold().eq("true").fillna(False)
    matched = sorted(
        original
        & set(rna.index)
        & set(pro.index)
        & set(cnv.index)
        & set(covered[covered].index)
    )
    if len(matched) < 400:
        raise RuntimeError(f"Unexpected eligible cohort size {len(matched)}")
    mats = {
        "RNA": rna.reindex(matched).astype("float32"),
        "PRO": pro.reindex(matched).astype("float32"),
        "CNV": cnv.reindex(matched).astype("float32"),
        "MUT": mut.reindex(matched).fillna(0).astype("float32"),
    }
    meta = splits.set_index("model_id").loc[matched]
    raw = pd.read_excel(
        h.datafile("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    raw["model_id"] = norm(raw.SANGER_MODEL_ID)
    raw["drug_id"] = norm(raw.DRUG_ID)
    raw.AUC = pd.to_numeric(raw.AUC, errors="coerce")
    raw = raw[
        raw.model_id.isin(matched) & raw.drug_id.isin(DRUGS) & raw.AUC.between(0, 1)
    ]
    if raw.duplicated(["drug_id", "model_id"]).any():
        raise ValueError("Duplicate compound-cellline records")
    say(
        "COHORT",
        f"Original TRAIN matched multi-omics models={len(matched)}. Only TRAIN responses used.",
    )
    # If a legacy similarity table is present, report its differing definition for
    # audit purposes. The current analysis does not require or reuse that table.
    a = ROOT / "Results" / "Phase2E" / "cancer_pairwise_response_correlations.csv"
    b = ROOT / "Results" / "Phase6A" / "focus_pairwise_drug_response_correlations.csv"
    for key, path in [("Phase2E", a), ("Phase6A", b)]:
        if path.exists():
            df = pd.read_csv(path)
            x = df[
                (
                    (df.cancer_a.astype(str).str.contains("Non-Small Cell Lung|NSCLC"))
                    & (df.cancer_b.astype(str).str.contains("Pancreatic"))
                )
                | (
                    (df.cancer_b.astype(str).str.contains("Non-Small Cell Lung|NSCLC"))
                    & (df.cancer_a.astype(str).str.contains("Pancreatic"))
                )
            ]
            if len(x):
                say(
                    "SIMILARITY AUDIT",
                    f"{key} pancreatic-vs-NSCLC reported rho={x.rho.iloc[0]:.4f}; n={x.iloc[0].get('shared_drugs', x.iloc[0].get('n_shared_drugs', '?'))}; different centering reference populations (all cancer categories vs focus cancers) and eligibility may change estimates.",
                )
    results = []
    oof = []
    checks = []
    for drug_id, name in DRUGS.items():
        d = raw[raw.drug_id == drug_id].set_index("model_id")
        ids = sorted(set(d.index) & set(matched))
        if len(ids) < 300:
            raise RuntimeError(f"{name}: only {len(ids)} eligible lines")
        loc = pd.Index(matched).get_indexer(ids)
        x = {k: v.iloc[loc].reset_index(drop=True) for k, v in mats.items()}
        y = d.loc[ids, "AUC"].to_numpy(dtype=float)
        can = meta.loc[ids, "cancer_type"].astype(str).to_numpy()
        groups = meta.loc[ids, "group_key"].astype(str).to_numpy()
        say(
            "DRUG",
            f"{name} [{drug_id}] training cohort={len(ids)}, colorectal={(can == TARGETS[0]).sum()}, pancreatic={(can == TARGETS[1]).sum()}",
        )
        unique = np.unique(groups)
        for repeat in range(1, REPEATS + 1):
            rng = np.random.default_rng(SEED + repeat)
            order = rng.permutation(len(unique))
            groupmap = dict(zip(unique, order))
            scrambled = np.array([groupmap[g] for g in groups])
            for fold, (ti, vi) in enumerate(
                GroupKFold(n_splits=FOLDS).split(ids, groups=scrambled), 1
            ):
                if set(groups[ti]) & set(groups[vi]):
                    raise AssertionError("Patient group leaked")
                selected = h.feature_sets(x, ti)
                # Same exact fold-wise molecular inputs for all three model types.
                blocks = ["RNA", "PRO", "MUT", "CNV"]
                tr = np.hstack(
                    [x[k].iloc[ti][selected[k]].to_numpy(np.float32) for k in blocks]
                )
                va = np.hstack(
                    [x[k].iloc[vi][selected[k]].to_numpy(np.float32) for k in blocks]
                )
                imputer = SimpleImputer(strategy="median", keep_empty_features=True)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    tr = imputer.fit_transform(tr).astype("float32")
                    va = imputer.transform(va).astype("float32")
                n_rna = len(selected["RNA"])
                n_pro = len(selected["PRO"])
                # Mean neighbors of RNA genes and of observed protein abundances separately.
                wm = []
                randm = []
                cover = []
                for k in ["RNA", "PRO"]:
                    W, nn, ee = network_matrix(selected[k], ppi)
                    Wr, _, _ = network_matrix(
                        selected[k], ppi, seed=SEED + repeat * 100 + fold
                    )
                    wm.append(W)
                    randm.append(Wr)
                    cover.append((k, nn, len(selected[k]), ee))
                if sum(z[3] for z in cover) < 20:
                    raise RuntimeError(
                        f"{name} network coverage too low for reliable comparison: {cover}"
                    )

                # Neighbour aggregation in original feature units, common fold-only imputation.
                def extended(X, matrices):
                    r = X[:, :n_rna]
                    p = X[:, n_rna : n_rna + n_pro]
                    return np.hstack([X, r @ matrices[0].T, p @ matrices[1].T]).astype(
                        "float32"
                    )

                variants = {
                    "AllFour": (tr, va),
                    "PPI": (extended(tr, wm), extended(va, wm)),
                    "ShuffledPPI": (extended(tr, randm), extended(va, randm)),
                }
                for label, (xt, xv) in variants.items():
                    fit = model(SEED + repeat * 100 + fold)
                    fit.fit(xt, y[ti])
                    pred = fit.predict(xv)
                    results.append(
                        dict(
                            drug_id=drug_id,
                            drug_name=name,
                            repeat=repeat,
                            fold=fold,
                            model=label,
                            n_train=len(ti),
                            n_valid=len(vi),
                            MAE=float(mean_absolute_error(y[vi], pred)),
                        )
                    )
                    for ii, actual, predicted in zip(vi, y[vi], pred):
                        oof.append(
                            dict(
                                drug_id=drug_id,
                                drug_name=name,
                                repeat=repeat,
                                fold=fold,
                                model=label,
                                model_id=ids[ii],
                                cancer_type=can[ii],
                                auc=float(actual),
                                pred=float(predicted),
                            )
                        )
                checks.append(
                    dict(
                        drug_id=drug_id,
                        repeat=repeat,
                        fold=fold,
                        rna_network_nodes=cover[0][1],
                        rna_selected=cover[0][2],
                        rna_edges=cover[0][3],
                        protein_network_nodes=cover[1][1],
                        protein_selected=cover[1][2],
                        protein_edges=cover[1][3],
                    )
                )
    scores = pd.DataFrame(results)
    pred = pd.DataFrame(oof)
    cov = pd.DataFrame(checks)
    scores.to_csv(OUT / "ppi_vs_multiomics_fold_scores.csv", index=False)
    pred.to_csv(OUT / "ppi_vs_multiomics_oof_predictions.csv", index=False)
    cov.to_csv(OUT / "ppi_fold_network_coverage.csv", index=False)
    summ = scores.groupby(["drug_id", "drug_name", "model"], as_index=False).MAE.mean()
    wide = scores.pivot(
        index=["drug_id", "drug_name", "repeat", "fold"], columns="model", values="MAE"
    ).reset_index()
    wide["gain_PPI_vs_baseline"] = wide.AllFour - wide.PPI
    wide["gain_PPI_vs_shuffled"] = wide.ShuffledPPI - wide.PPI
    wide.to_csv(OUT / "ppi_paired_fold_differences.csv", index=False)
    per = wide.groupby(["drug_name"], as_index=False).agg(
        gain_vs_baseline=("gain_PPI_vs_baseline", "mean"),
        gain_vs_shuffle=("gain_PPI_vs_shuffled", "mean"),
        folds_better=("gain_PPI_vs_baseline", lambda x: int((x > 0).sum())),
    )
    sub = (
        pred[pred.cancer_type.isin(TARGETS)]
        .groupby(["drug_name", "cancer_type", "repeat", "model"])
        .apply(
            lambda t: pd.Series(
                {"n": len(t), "MAE": mean_absolute_error(t.auc, t.pred)}
            ),
            include_groups=False,
        )
        .reset_index()
    )
    sub.to_csv(OUT / "ppi_target_cancer_repeat_oof_MAE.csv", index=False)
    audit = (
        pred.groupby(["drug_name", "repeat", "model", "cancer_type"])
        .size()
        .unstack(fill_value=0)
    )
    audit.to_csv(OUT / "cancer_type_fold_coverage.csv")
    say(
        "NETWORK COVERAGE",
        "\n"
        + cov.groupby("drug_id")[
            ["rna_network_nodes", "rna_edges", "protein_network_nodes", "protein_edges"]
        ]
        .mean()
        .round(1)
        .to_string(),
    )
    say(
        "OVERALL SCORES",
        "\n"
        + summ.pivot(index="drug_name", columns="model", values="MAE")
        .round(5)
        .to_string(),
    )
    say("PPI PAIRED", "\n" + per.round(5).to_string(index=False))
    say(
        "TARGET CANCER OOF",
        "\n"
        + sub.groupby(["drug_name", "cancer_type", "model"])
        .agg(mean_MAE=("MAE", "mean"), n=("n", "mean"))
        .reset_index()
        .round(5)
        .to_string(index=False),
    )
    say(
        "CAUTION",
        "PPI-derived network features do not establish causal mechanisms; repeated CV folds overlap, comparisons are exploratory; colorectal/pancreatic subgroup sizes can be small.",
    )
    say(
        "CAUTION",
        "Shuffled node-label control disrupts gene-to-network assignment but does not preserve each gene identity’s degree. Network comparison does not demonstrate clinical efficacy or new indications.",
    )
    say("STATUS", f"SUCCESS elapsed={(time.time() - START) / 60:.1f} min; output={OUT}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        say("ERROR", f"{type(exc).__name__}: {exc}")
        say("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 95 + "\n=== BEGIN PROJECT 5 PHASE 6B RESULTS ===")
        print("\n".join(LOG))
        print("=== END PROJECT 5 PHASE 6B RESULTS ===\n" + "=" * 95)

# ----------------------------------------------------------------------------
# MAPK molecular association follow-up
# ----------------------------------------------------------------------------

from pathlib import Path
import importlib.util, warnings, time, traceback, zipfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase6C"
OUT.mkdir(parents=True, exist_ok=True)
LOG = []
START = time.time()
CANCERS = [
    "Colorectal Carcinoma",
    "Pancreatic Carcinoma",
    "Breast Carcinoma",
    "Non-Small Cell Lung Carcinoma",
    "Small Cell Lung Carcinoma",
]
DRUGS = {
    "1372": "Trametinib",
    "1060": "PD0325901",
    "1014": "Refametinib",
    "1564": "SCH772984",
    "1062": "Selumetinib",
    "1736": "Selumetinib (other ID)",
    "1908": "Ulixertinib",
}
# Prespecified core MAPK nodes, feedback regulators and context markers. Gene abundance != pathway activation.
GENES = [
    "KRAS",
    "NRAS",
    "HRAS",
    "BRAF",
    "RAF1",
    "ARAF",
    "MAP2K1",
    "MAP2K2",
    "MAPK1",
    "MAPK3",
    "DUSP4",
    "DUSP6",
    "SPRY2",
    "SPRY4",
    "ETV4",
    "ETV5",
    "FOS",
    "JUN",
    "MYC",
    "EGFR",
    "ERBB2",
    "NF1",
    "PTEN",
    "PIK3CA",
    "AKT1",
    "CDKN2A",
    "TP53",
    "GAPDH",
]
SIGNATURE = [
    "DUSP4",
    "DUSP6",
    "SPRY2",
    "SPRY4",
    "ETV4",
    "ETV5",
]  # transcriptional feedback proxy, not direct MAPK activity
MUT_GENES = ["KRAS", "NRAS", "BRAF", "NF1", "TP53", "PIK3CA", "APC", "SMAD4"]
MIN_N = 12


def say(tag, v):
    LOG.append(f"[{tag}] {v}")


def norm(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def file(prefix):
    a = list(DATA.glob(prefix + "*"))
    if len(a) != 1:
        raise FileNotFoundError(
            f"Expected one Data/{prefix}*; found {[p.name for p in a]}"
        )
    return a[0]


def disp(df, n=16):
    return (
        "No eligible results"
        if df.empty
        else df.head(n).to_string(index=False, float_format=lambda x: f"{x:.4g}")
    )


def helper():
    p = ROOT / "04_multiomics_explainability.py"
    if not p.exists():
        raise FileNotFoundError(
            "Phase3A script needed to read existing molecular tables"
        )
    spec = importlib.util.spec_from_file_location("phase3a_mapk_helper", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def read_rna(models):
    path = file("rnaseq_merged_")
    with zipfile.ZipFile(path) as z:
        names = [
            n for n in z.namelist() if "rsem_tpm_" in n.lower() and n.endswith(".csv")
        ]
        if len(names) != 1:
            raise ValueError(f"Unexpected RNA archive contents: {names}")
        with z.open(names[0]) as fh:
            hdr = pd.read_csv(fh, nrows=0).columns.tolist()
        cols = [c for c in hdr if c in models and c.startswith("SIDM")]
        if len(cols) < 100:
            raise ValueError(f"RNA model columns: {len(cols)}")
        with z.open(names[0]) as fh:
            tmp = pd.read_csv(fh, usecols=[hdr[0]] + cols, low_memory=False)
    if tmp.iloc[:3, 0].astype(str).tolist() != [
        "model_name",
        "data_source",
        "gene_symbol",
    ]:
        raise ValueError("RNA metadata-row structure changed")
    sources = tmp.loc[1, cols].astype(str).rename("rna_source")
    part = tmp.iloc[3:].copy()
    part["symbol"] = part.iloc[:, 0].astype(str).str.upper().str.strip()
    part = part[part.symbol.isin(GENES)]
    if part.symbol.duplicated().any():
        say(
            "RNA",
            f"Excluded duplicated symbols: {part.loc[part.symbol.duplicated(False), 'symbol'].unique().tolist()}",
        )
        part = part[~part.symbol.duplicated(False)]
    if part.empty:
        raise ValueError("MAPK genes not found in RNA matrix")
    arr = part[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype="float32")
    arr[(arr < 0) | ~np.isfinite(arr)] = np.nan
    arr = np.log2(arr + 1)
    rna = pd.DataFrame(arr.T, index=cols, columns=part.symbol.tolist())
    rna.index.name = "model_id"
    rna["rna_source"] = rna.index.map(sources.to_dict())
    say(
        "RNA",
        f"Matched models={rna.shape[0]}; selected genes={rna.shape[1] - 1}; sources={rna.rna_source.value_counts(dropna=False).to_dict()}",
    )
    return rna


def stats_assoc(df, drug, cancer, modality, feature, typ):
    d = df[["AUC", feature]].dropna()
    n = len(d)
    if typ == "continuous":
        if n < MIN_N or d[feature].nunique() < 4 or d.AUC.nunique() < 4:
            return None
        r, p = spearmanr(d[feature], d.AUC)
        if not np.isfinite(p):
            return None
        return dict(
            drug=drug,
            cancer=cancer,
            modality=modality,
            feature=feature,
            n=n,
            event_n=np.nan,
            effect=float(r),
            p=float(p),
            metric="Spearman rho",
        )
    if typ == "mutation":
        if n < MIN_N:
            return None
        a = d.loc[d[feature].eq(1), "AUC"]
        b = d.loc[d[feature].eq(0), "AUC"]
        if min(len(a), len(b)) < 4:
            return None
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        return dict(
            drug=drug,
            cancer=cancer,
            modality=modality,
            feature=feature,
            n=n,
            event_n=len(a),
            effect=float(a.median() - b.median()),
            p=float(p),
            metric="median AUC event-minus-no-record",
        )


def rank_adjusted(df, feature):
    # Joint colorectal/pancreatic descriptive regression: conditional association, HC3 robust SE.
    d = (
        df[["AUC", feature, "cancer_type", "rna_source"]]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )
    if len(d) < 35 or d[feature].nunique() < 8 or d.cancer_type.nunique() < 2:
        return None
    d = d.copy()
    d["rx"] = d[feature].rank(pct=True)
    d["ry"] = d.AUC.rank(pct=True)
    d["is_panc"] = (d.cancer_type == "Pancreatic Carcinoma").astype(float)
    d["interaction"] = d.rx * d.is_panc
    X = pd.concat(
        [
            d[["rx", "is_panc", "interaction"]],
            pd.get_dummies(d.rna_source, prefix="src", drop_first=True, dtype=float),
        ],
        axis=1,
    )
    X = sm.add_constant(X, has_constant="add")
    try:
        res = sm.OLS(d.ry, X.astype(float)).fit(cov_type="HC3")
    except (ValueError, np.linalg.LinAlgError):
        return None
    return dict(
        feature=feature,
        n=len(d),
        n_crc=int((d.is_panc == 0).sum()),
        n_panc=int((d.is_panc == 1).sum()),
        crc_rank_slope=float(res.params["rx"]),
        p_crc_slope=float(res.pvalues["rx"]),
        pancreatic_delta_slope=float(res.params["interaction"]),
        p_interaction=float(res.pvalues["interaction"]),
        pancreatic_rank_slope=float(res.params["rx"] + res.params["interaction"]),
    )


def main():
    say(
        "DESIGN",
        "Observational GDSC2 MAPK molecular associations in colorectal and pancreatic cancers; breast/NSCLC/SCLC comparator groups.",
    )
    say(
        "DESIGN",
        "Drug-specific AUC comparisons only. Mutation summary event/no-record ≠ functional mutant/wildtype; RNA/protein abundance ≠ MAPK activity.",
    )
    say(
        "DESIGN",
        "Within-cancer associations and cancer/source-adjusted joint models are exploratory; BH-FDR across all reported tests in each family.",
    )
    models = pd.read_csv(file("model_list_"), low_memory=False)
    say("MODEL SCHEMA", f"Model metadata cols={list(models.columns)[:18]}")
    raw = pd.read_excel(
        file("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    raw["model_id"] = norm(raw.SANGER_MODEL_ID)
    raw["drug_id"] = norm(raw.DRUG_ID)
    raw.AUC = pd.to_numeric(raw.AUC, errors="coerce")
    raw = raw[raw.drug_id.isin(DRUGS) & raw.AUC.between(0, 1)].copy()
    # Use existing original Phase2 cancer annotation; no reclassification or silent broad grouping.
    splits = pd.read_csv(
        ROOT / "Results" / "Phase2" / "model_splits.csv", dtype={"model_id": "string"}
    )
    splits.model_id = norm(splits.model_id)
    if splits.model_id.duplicated().any():
        raise ValueError("Nonunique original model identifiers")
    meta = splits[["model_id", "cancer_type", "split"]].copy()
    df = raw.merge(meta, on="model_id", how="inner", validate="many_to_one")
    df = df[df.cancer_type.isin(CANCERS)]
    if df.duplicated(["model_id", "drug_id"]).any():
        raise ValueError("Duplicate model-drug observations")
    say(
        "COHORT",
        f"{df.model_id.nunique()} matched target/comparator models; {len(df)} drug observations; 7 drug IDs; counts="
        + str(df.groupby("cancer_type").model_id.nunique().to_dict()),
    )
    rna = read_rna(set(df.model_id))
    p3 = helper()
    prot = p3.read_proteome(set(df.model_id))
    cnv = p3.read_cnv(set(df.model_id))
    mut = p3.read_mutation(set(df.model_id))
    # Preserve NA in proteomics/CNV; mutation summary non-observation is marked zero ONLY for WES coverage-flagged models.
    cov = pd.read_csv(file("model_dataset_availability_"), low_memory=False)
    cov.model_id = norm(cov.model_id)
    flags = [
        c
        for c in ("Mutation Sanger Cell Lines (WES)", "Mutation Broad Cell Lines (WES)")
        if c in cov
    ]
    if not flags:
        raise ValueError("Mutation WES coverage flags unavailable")
    covered = pd.Series(False, index=cov.index)
    for c in flags:
        covered |= cov[c].astype("string").str.lower().eq("true").fillna(False)
    cover = set(cov.loc[covered, "model_id"])
    say(
        "MODALITIES",
        f"RNA={len(rna)}, proteome={len(prot)}, CNV={len(cnv)}, WES coverage={len(set(df.model_id) & cover)}",
    )
    # Add selected expression, proteome and CNV variables with explicit prefixes.
    expression = [g for g in GENES if g in rna]
    protein = [g for g in GENES if g in prot]
    copies = [g for g in MUT_GENES if g in cnv]
    merged = df.merge(
        rna[["rna_source"] + expression]
        .rename(columns={g: "RNA:" + g for g in expression})
        .reset_index(),
        on="model_id",
        how="left",
        validate="many_to_one",
    )
    pro = prot.reindex(columns=protein).rename(columns={g: "PRO:" + g for g in protein})
    pro.index.name = "model_id"
    merged = merged.merge(
        pro.reset_index(), on="model_id", how="left", validate="many_to_one"
    )
    cop = cnv.reindex(columns=copies).rename(columns={g: "CNV:" + g for g in copies})
    cop.index.name = "model_id"
    merged = merged.merge(
        cop.reset_index(), on="model_id", how="left", validate="many_to_one"
    )
    mutation = pd.DataFrame(index=sorted(set(merged.model_id) & cover))
    for g in MUT_GENES:
        mutation["MUT:" + g] = (
            mut.reindex(index=mutation.index, columns=[g]).fillna(0)[g].astype(float)
        )
    mutation.index.name = "model_id"
    merged = merged.merge(
        mutation.reset_index(), on="model_id", how="left", validate="many_to_one"
    )
    # Within-source standardized transcriptional feedback score; only score if >=4/6 measured.
    parts = []
    for src, d in merged.groupby("rna_source", dropna=True):
        x = d[["RNA:" + g for g in SIGNATURE if "RNA:" + g in merged]].copy()
        if x.shape[1] < 4:
            continue
        mean = x.mean()
        std = x.std(ddof=1).replace(0, np.nan)
        z = (x - mean) / std
        score = z.mean(axis=1).where(z.notna().sum(axis=1) >= 4)
        parts.append(score.rename("RNA:MAPK_feedback_score"))
    merged["RNA:MAPK_feedback_score"] = (
        pd.concat(parts).reindex(merged.index) if parts else np.nan
    )
    say(
        "FEATURES",
        f"RNA genes={len(expression)}, protein genes={len(protein)}, CNV genes={len(copies)}, mutation summary genes={len(MUT_GENES)}; feedback-score genes={SIGNATURE}",
    )
    merged[
        ["model_id", "drug_id", "DRUG_NAME", "cancer_type", "AUC", "rna_source"]
        + ["RNA:" + g for g in expression]
        + ["PRO:" + g for g in protein]
        + ["CNV:" + g for g in copies]
        + ["MUT:" + g for g in MUT_GENES]
        + ["RNA:MAPK_feedback_score"]
    ].to_csv(
        OUT / "target_drug_cancer_molecular_matrix.csv.gz",
        index=False,
        compression="gzip",
    )
    rows = []
    for (drug, cancer), d in merged.groupby(["drug_id", "cancer_type"]):
        for m, columns, typ in [
            (
                "RNA",
                ["RNA:" + g for g in expression] + ["RNA:MAPK_feedback_score"],
                "continuous",
            ),
            ("PRO", ["PRO:" + g for g in protein], "continuous"),
            ("CNV", ["CNV:" + g for g in copies], "continuous"),
            ("MUT", ["MUT:" + g for g in MUT_GENES], "mutation"),
        ]:
            for g in columns:
                r = stats_assoc(d, DRUGS[drug], cancer, m, g, typ)
                if r is not None:
                    r.update(drug_id=drug)
                    rows.append(r)
    a = pd.DataFrame(rows)
    if a.empty:
        raise RuntimeError("No eligible molecular association tests")
    a["q_fdr"] = multipletests(a.p, method="fdr_bh")[1]
    a = a.sort_values(["q_fdr", "p"])
    a.to_csv(OUT / "all_within_cancer_molecular_associations.csv", index=False)
    say(
        "TESTS",
        f"Within-cancer association tests={len(a)}; q<.05={(a.q_fdr < 0.05).sum()}; target-cancer significant={int(((a.q_fdr < 0.05) & a.cancer.isin(CANCERS[:2])).sum())}. Mutation classes cannot be inferred from summary.",
    )
    say(
        "TARGET TOP",
        "\n"
        + disp(
            a[a.cancer.isin(CANCERS[:2])][
                [
                    "drug",
                    "cancer",
                    "modality",
                    "feature",
                    "n",
                    "event_n",
                    "effect",
                    "q_fdr",
                ]
            ],
            25,
        ),
    )
    say(
        "MAPK SCORE",
        "\n"
        + disp(
            a[a.feature.eq("RNA:MAPK_feedback_score")][
                ["drug", "cancer", "n", "effect", "p", "q_fdr"]
            ],
            20,
        ),
    )
    joint = []
    for drug, d in merged[merged.cancer_type.isin(CANCERS[:2])].groupby("drug_id"):
        for feature in ["RNA:" + g for g in GENES if g in expression] + [
            "RNA:MAPK_feedback_score"
        ]:
            z = rank_adjusted(d, feature)
            if z is not None:
                z.update(drug=DRUGS[drug], drug_id=drug)
                joint.append(z)
    b = pd.DataFrame(joint)
    if not b.empty:
        b["q_interaction"] = multipletests(b.p_interaction, method="fdr_bh")[1]
        b["q_crc_slope"] = multipletests(b.p_crc_slope, method="fdr_bh")[1]
        b.sort_values("q_interaction").to_csv(
            OUT / "crc_pancreatic_rna_interaction_tests.csv", index=False
        )
        say(
            "INTERACTIONS",
            f"Tested {len(b)} cancer×expression interactions; FDR<.05={(b.q_interaction < 0.05).sum()}. HC3 robust association, not a causal effect.",
        )
        say(
            "INTERACTION TOP",
            "\n"
            + disp(
                b.sort_values("q_interaction")[
                    [
                        "drug",
                        "feature",
                        "n_crc",
                        "n_panc",
                        "crc_rank_slope",
                        "pancreatic_rank_slope",
                        "pancreatic_delta_slope",
                        "q_interaction",
                    ]
                ],
                12,
            ),
        )
    else:
        say("INTERACTIONS", "No eligible joint RNA analyses")
    summary = (
        merged.groupby(["drug_id", "DRUG_NAME", "cancer_type"])
        .agg(
            n=("AUC", "size"),
            median_auc=("AUC", "median"),
            iqr_auc=("AUC", lambda s: s.quantile(0.75) - s.quantile(0.25)),
        )
        .reset_index()
    )
    summary.to_csv(OUT / "drug_cancer_response_summary.csv", index=False)
    say(
        "RESPONSE TARGETS",
        "\n" + disp(summary[summary.cancer_type.isin(CANCERS[:2])], 22),
    )
    say(
        "INTERPRETATION",
        "Positive Spearman rho: higher abundance is associated with higher AUC (usually weaker experimental viability suppression); negative is the opposite. Mutation effect is event-recorded minus no-record median AUC.",
    )
    say(
        "LIMITS",
        "27 pancreatic models in Phase6A, fewer in omics-complete subsets; correlated predictors and multiple treatments share cell lines; cell lines are not patients. Screening AUC is not clinical efficacy.",
    )
    say(
        "LIMITS",
        "Source-adjusted rank interaction fits are exploratory. Source effects, KRAS variant identity, mutation record completeness, protein missingness and pathway-feedback score validity require caution. No claims of new indications.",
    )
    say(
        "STATUS",
        f"SUCCESS elapsed={(time.time() - START) / 60:.1f} min; saved to {OUT}",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        say("ERROR", f"{type(e).__name__}: {e}")
        say("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 98 + "\n=== BEGIN PROJECT 5 PHASE 6C RESULTS ===")
        print("\n".join(LOG))
        print("=== END PROJECT 5 PHASE 6C RESULTS ===\n" + "=" * 98)

# ----------------------------------------------------------------------------
# MAPK signature incremental prediction
# ----------------------------------------------------------------------------

from pathlib import Path
import traceback, time, warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score

ROOT = PROJECT_ROOT
OUT = ROOT / "Results" / "Phase6D"
OUT.mkdir(parents=True, exist_ok=True)
INFILE = ROOT / "Results" / "Phase6C" / "target_drug_cancer_molecular_matrix.csv.gz"
SPLITFILE = ROOT / "Results" / "Phase2" / "model_splits.csv"
CANCERS = ["Colorectal Carcinoma", "Pancreatic Carcinoma"]
DRUGS = {
    "1372": "Trametinib",
    "1060": "PD0325901",
    "1014": "Refametinib",
    "1564": "SCH772984",
    "1062": "Selumetinib (1062)",
    "1736": "Selumetinib (1736)",
    "1908": "Ulixertinib",
}
SIG = ["RNA:" + g for g in ["DUSP4", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5"]]
MUT = ["MUT:" + g for g in ["KRAS", "BRAF", "NRAS", "NF1"]]
MODELS = {
    "Cancer+Source": (False, False),
    "Cancer+Source+Mut": (True, False),
    "Cancer+Source+Signature": (False, True),
    "Cancer+Source+Both": (True, True),
}
SEED = 20261028
REPEATS = 3
FOLDS = 4
LOG = []
T0 = time.time()


def say(k, x):
    LOG.append(f"[{k}] {x}")


def show(df, n=30):
    return (
        "No eligible results"
        if df.empty
        else df.head(n).to_string(index=False, float_format=lambda x: f"{x:.5g}")
    )


def normal(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def signature_fold(train, val):
    # Each gene normalised using training data, score is mean z-value over >=4 genes.
    a = train[SIG].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    b = val[SIG].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    means = np.nanmean(a, axis=0)
    sds = np.nanstd(a, axis=0, ddof=1)
    valid = np.isfinite(means) & np.isfinite(sds) & (sds > 1e-8)
    if valid.sum() < 4:
        raise ValueError("Fewer than 4 variable feedback genes in training fold")

    def score(x):
        z = (x[:, valid] - means[valid]) / sds[valid]
        num = np.isfinite(z).sum(axis=1)
        return np.where(num >= 4, np.nansum(z, axis=1) / np.maximum(num, 1), np.nan)

    return score(a), score(b), int(valid.sum())


def matrices(tr, va, with_mut, with_sig, st, sv):
    # Same nuisance covariates for ALL models; cancer-aware reference is not train mean.
    def base(d):
        return np.column_stack(
            [
                d.cancer_type.eq("Pancreatic Carcinoma").to_numpy(dtype=float),
                d.rna_source.astype("string")
                .fillna("Missing")
                .eq("Broad")
                .to_numpy(dtype=float),
                d.rna_source.astype("string")
                .fillna("Missing")
                .eq("Sanger")
                .to_numpy(dtype=float),
            ]
        )

    x = base(tr)
    y = base(va)
    if with_mut:
        x = np.column_stack(
            [x, tr[MUT].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)]
        )
        y = np.column_stack(
            [y, va[MUT].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)]
        )
    if with_sig:
        x = np.column_stack([x, st])
        y = np.column_stack([y, sv])
    return x, y


def run():
    say(
        "DESIGN",
        "Original Phase2 TRAIN cell lines only; 3x4 stratified patient-group CV, identical held-out samples per model.",
    )
    say(
        "DESIGN",
        "Fixed Ridge(alpha=10) for 7 selected MAPK drug IDs; score z-normalisation, imputation, scaling fitted within each fold.",
    )
    say(
        "DESIGN",
        "All models adjust for cancer label and RNA data source; mutation records mean event-present vs no recorded event, not mutant/wildtype.",
    )
    if not INFILE.exists():
        raise FileNotFoundError(f"Run Phase6C first; missing {INFILE}")
    if not SPLITFILE.exists():
        raise FileNotFoundError(SPLITFILE)
    d = pd.read_csv(
        INFILE, low_memory=False, dtype={"drug_id": "string", "model_id": "string"}
    )
    s = pd.read_csv(SPLITFILE, low_memory=False, dtype={"model_id": "string"})
    for col in ("group_key", "split"):
        if col not in s:
            raise ValueError(f"Missing Phase2 field: {col}")
    d.model_id = normal(d.model_id)
    s.model_id = normal(s.model_id)
    d.drug_id = normal(d.drug_id)
    if s.model_id.duplicated().any():
        raise ValueError("Nonunique model IDs in Phase2 splits")
    required = ["model_id", "drug_id", "cancer_type", "AUC", "rna_source"] + SIG + MUT
    missing = [c for c in required if c not in d]
    if missing:
        raise ValueError(f"Missing molecular matrix columns: {missing}")
    d = d.merge(
        s[["model_id", "group_key", "split"]],
        on="model_id",
        how="inner",
        validate="many_to_one",
    )
    d = d[
        d.split.eq("train") & d.cancer_type.isin(CANCERS) & d.drug_id.isin(DRUGS)
    ].copy()
    d["AUC"] = pd.to_numeric(d.AUC, errors="coerce")
    d = d[d.AUC.between(0, 1)]
    if d.duplicated(["model_id", "drug_id"]).any():
        raise ValueError("Duplicate drug-cell-line observations")
    if d.group_key.isna().any():
        raise ValueError("Missing patient grouping")
    say(
        "DATA",
        f"{d.model_id.nunique()} unique original TRAIN cell lines; {len(d)} drug observations; "
        + str(d.groupby("cancer_type").model_id.nunique().to_dict()),
    )
    rows = []
    predictions = []
    feature_counts = []
    for did, name in DRUGS.items():
        drug = d[d.drug_id.eq(did)].copy().reset_index(drop=True)
        counts = drug.groupby("cancer_type").model_id.nunique().to_dict()
        if min(counts.get(c, 0) for c in CANCERS) < 12:
            say("SKIP", f"{name} ID={did} insufficient group sample sizes: {counts}")
            continue
        say("DRUG", f"{name} [{did}]: n={len(drug)}, cohorts={counts}")
        for rep in range(REPEATS):
            cv = StratifiedGroupKFold(
                n_splits=FOLDS, shuffle=True, random_state=SEED + rep
            )
            for fold, (ti, vi) in enumerate(
                cv.split(drug, drug.cancer_type, drug.group_key), 1
            ):
                tr = drug.iloc[ti]
                va = drug.iloc[vi]
                if set(tr.group_key) & set(va.group_key):
                    raise ValueError("Group leakage")
                if tr.cancer_type.nunique() != 2 or va.cancer_type.nunique() != 2:
                    raise ValueError(
                        f"Nonrepresentative fold for {name}, repeat {rep}, fold {fold}; consider fewer folds"
                    )
                st, sv, n_genes = signature_fold(tr, va)
                feature_counts.append(
                    dict(drug_id=did, repeat=rep + 1, fold=fold, score_genes=n_genes)
                )
                yy = tr.AUC.to_numpy(float)
                true = va.AUC.to_numpy(float)
                for label, (addmut, addsig) in MODELS.items():
                    a, b = matrices(tr, va, addmut, addsig, st, sv)
                    imp = SimpleImputer(strategy="median", keep_empty_features=True)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        a = imp.fit_transform(a)
                        b = imp.transform(b)
                    sc = StandardScaler()
                    a = sc.fit_transform(a)
                    b = sc.transform(b)
                    model = Ridge(alpha=10.0)
                    model.fit(a, yy)
                    predicted = np.clip(model.predict(b), 0, 1)
                    rows.append(
                        dict(
                            drug_id=did,
                            drug=name,
                            repeat=rep + 1,
                            fold=fold,
                            model=label,
                            n_train=len(tr),
                            n_val=len(va),
                            MAE=mean_absolute_error(true, predicted),
                        )
                    )
                    for ix, (obs, est) in enumerate(zip(true, predicted)):
                        predictions.append(
                            dict(
                                drug_id=did,
                                drug=name,
                                repeat=rep + 1,
                                fold=fold,
                                model=label,
                                model_id=va.model_id.iloc[ix],
                                cancer_type=va.cancer_type.iloc[ix],
                                observed_AUC=obs,
                                predicted_AUC=est,
                            )
                        )
    folds = pd.DataFrame(rows)
    pred = pd.DataFrame(predictions)
    if folds.empty:
        raise RuntimeError("No valid drug CV evaluations")
    folds.to_csv(OUT / "fold_model_comparison.csv", index=False)
    pred.to_csv(OUT / "out_of_fold_predictions.csv.gz", index=False, compression="gzip")
    pd.DataFrame(feature_counts).to_csv(
        OUT / "fold_signature_gene_counts.csv", index=False
    )
    perf = folds.groupby(["drug_id", "drug", "model"], as_index=False).agg(
        mean_MAE=("MAE", "mean"), sd_MAE=("MAE", "std")
    )
    wide = perf.pivot(
        index=["drug_id", "drug"], columns="model", values="mean_MAE"
    ).reset_index()
    base = "Cancer+Source"
    both = "Cancer+Source+Both"
    mutate = "Cancer+Source+Mut"
    sig = "Cancer+Source+Signature"
    wide["gain_signature_vs_reference"] = wide[base] - wide[sig]
    wide["gain_mutation_vs_reference"] = wide[base] - wide[mutate]
    wide["gain_signature_beyond_mutation"] = wide[mutate] - wide[both]
    wide["gain_both_vs_reference"] = wide[base] - wide[both]
    wide.to_csv(OUT / "per_drug_incremental_gains.csv", index=False)
    bytype = (
        pred.groupby(["drug_id", "drug", "cancer_type", "model"])
        .apply(
            lambda x: pd.Series(
                {
                    "n_unique_lines": x.model_id.nunique(),
                    "MAE": np.mean(np.abs(x.observed_AUC - x.predicted_AUC)),
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )
    bytype.to_csv(OUT / "per_cancer_out_of_fold_MAE.csv", index=False)
    ct = bytype.pivot(
        index=["drug_id", "drug", "cancer_type", "n_unique_lines"],
        columns="model",
        values="MAE",
    ).reset_index()
    ct["gain_signature_beyond_mutation"] = ct[mutate] - ct[both]
    ct.to_csv(OUT / "per_cancer_incremental_gains.csv", index=False)
    # Descriptive only: CV repeats reuse cell lines, not independent samples.
    say(
        "MODEL SUMMARY",
        "\n"
        + show(
            perf.groupby("model", as_index=False).agg(
                mean_drug_MAE=("mean_MAE", "mean"),
                median_drug_MAE=("mean_MAE", "median"),
            )
        ),
    )
    say(
        "PER DRUG",
        "\n"
        + show(
            wide[
                [
                    "drug",
                    "gain_signature_vs_reference",
                    "gain_mutation_vs_reference",
                    "gain_signature_beyond_mutation",
                    "gain_both_vs_reference",
                ]
            ],
            20,
        ),
    )
    say(
        "PER CANCER",
        "\n"
        + show(
            ct[
                [
                    "drug",
                    "cancer_type",
                    "n_unique_lines",
                    "gain_signature_beyond_mutation",
                ]
            ],
            20,
        ),
    )
    say(
        "SIGNATURE",
        f"Feedback score gene count per fold: {pd.DataFrame(feature_counts).score_genes.value_counts().to_dict()}",
    )
    say(
        "SUMMARY",
        f"Signature improves over mutation-adjusted model in {(wide.gain_signature_beyond_mutation > 0).sum()}/{len(wide)} drugs; descriptive mean paired gain={wide.gain_signature_beyond_mutation.mean():+.5f}",
    )
    say(
        "LIMITATIONS",
        "Only approximately 30 CRC and 19 pancreatic original TRAIN lines per drug. Repeated folds overlap; no independent confirmation or clinical prediction.",
    )
    say(
        "LIMITATIONS",
        "Recorded mutation events do not establish mutation function; MAPK feedback score is RNA proxy, not phospho-ERK activity. Selected compounds and signature were motivated by earlier observations.",
    )
    say(
        "CAUTION",
        "Models adjust for RNA source with a binary nuisance term; small cohorts may produce unstable estimates. Same selected drugs and source cohort; no PRISM outcome evaluation.",
    )
    say("STATUS", f"SUCCESS elapsed={time.time() - T0:.1f}s; outputs={OUT}")


if __name__ == "__main__":
    try:
        run()
    except Exception as exc:
        say("ERROR", f"{type(exc).__name__}: {exc}")
        say("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 100 + "\n=== BEGIN PROJECT 5 PHASE 6D RESULTS ===")
        print("\n".join(LOG))
        print("=== END PROJECT 5 PHASE 6D RESULTS ===\n" + "=" * 100)
