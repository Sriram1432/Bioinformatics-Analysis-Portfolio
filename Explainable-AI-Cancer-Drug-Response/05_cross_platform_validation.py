"""Audit GDSC2 and PRISM concordance and evaluate SLFN11 associations.

Compound matches remain assay-aware and adjusted rank analyses account for cancer lineage and RNA source where supported.

Run from the repository root with ``python 05_cross_platform_validation.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# Cross-platform compound and sample audit
# ----------------------------------------------------------------------------

from pathlib import Path
import time, traceback, re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase4A"
OUT.mkdir(parents=True, exist_ok=True)
REP = []
START = time.time()
FOCUS = {"1047": "Nutlin-3a (-)", "1003": "Camptothecin", "1060": "PD0325901"}


def say(k, s):
    REP.append(f"[{k}] {s}")


def file(fragment):
    a = list(DATA.glob(fragment + "*"))
    if len(a) != 1:
        raise FileNotFoundError(
            f"Expected one Data/{fragment}*; found {[x.name for x in a]}"
        )
    return a[0]


def col(df, candidates, required=True):
    m = {re.sub(r"[^a-z0-9]", "", str(k).lower()): k for k in df.columns}
    for c in candidates:
        k = m.get(re.sub(r"[^a-z0-9]", "", c.lower()))
        if k is not None:
            return k
    if required:
        raise KeyError(f"None of {candidates} in columns: {list(df.columns)}")
    return None


def ident(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def names(s):
    return (
        s.astype("string")
        .str.casefold()
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
        .str.replace("−", "-", regex=False)
    )


def pretty(df, cols, n=12):
    return (
        "No eligible records"
        if df.empty
        else df[cols].head(n).to_string(index=False, float_format=lambda x: f"{x:.5g}")
    )


def main():
    say(
        "DESIGN",
        "Cross-platform assay concordance on EXACT matched cell lines; GDSC2 and PRISM AUC may have different scales/dose ranges.",
    )
    say(
        "DESIGN",
        "Name-only matched drugs are PROVISIONAL, not chemical identity verification. No confirmatory efficacy, repurposing or external model claims.",
    )
    g = pd.read_excel(
        file("GDSC2_fitted_dose_response_"),
        usecols=["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"],
    )
    g["model_id"] = ident(g.SANGER_MODEL_ID)
    g["drug_id"] = ident(g.DRUG_ID)
    g["auc_gdsc"] = pd.to_numeric(g.AUC, errors="coerce")
    g = g[g.auc_gdsc.between(0, 1) & g.model_id.notna()].copy()
    if g.duplicated(["drug_id", "model_id"]).any():
        raise ValueError("GDSC has duplicate drug_id × model_id rows")
    model = pd.read_csv(file("model_list_"), low_memory=False)
    mc = col(model, ["model_id"])
    bc = col(model, ["BROAD_ID", "DepMap_ID", "depmap_id", "broad_id"], required=False)
    if bc is None:
        raise KeyError(f"No Broad/DepMap ID field in model list: {list(model.columns)}")
    model["model_id"] = ident(model[mc])
    model["depmap_id"] = ident(model[bc])
    model = model[["model_id", "depmap_id"]].dropna().drop_duplicates()
    if (
        model.groupby("model_id").depmap_id.nunique().max() > 1
        or model.groupby("depmap_id").model_id.nunique().max() > 1
    ):
        raise ValueError("Ambiguous model ID to DepMap ID mapping; inspect model_list_")
    g = g.merge(model, on="model_id", how="left", validate="many_to_one")
    splpath = ROOT / "Results" / "Phase2" / "model_splits.csv"
    if splpath.exists():
        sp = pd.read_csv(splpath)
        sp["model_id"] = ident(sp.model_id)
        g = g.merge(
            sp[["model_id", "split"]], on="model_id", how="left", validate="many_to_one"
        )
    else:
        g["split"] = pd.NA
    pp = file("secondary-screen-dose-response-curve-parameters")
    hdr = pd.read_csv(pp, nrows=0).columns.tolist()
    say("SCHEMA", f"PRISM columns={hdr}")
    required = ["broad_id", "depmap_id", "name", "auc"]
    optional = [
        "r2",
        "screen_id",
        "smiles",
        "passed_str_profiling",
        "moa",
        "target",
        "indication",
    ]
    use = [c for c in required + optional if c in hdr]
    if not set(required).issubset(use):
        raise KeyError(
            "PRISM secondary dose response missing required columns; no guesses made"
        )
    chunks = []
    rows = 0
    for ch in pd.read_csv(pp, usecols=use, chunksize=150000, low_memory=False):
        rows += len(ch)
        ch["depmap_id"] = ident(ch.depmap_id)
        ch["broad_id"] = ident(ch.broad_id)
        ch["auc_prism"] = pd.to_numeric(ch.auc, errors="coerce")
        ch = ch[
            ch.auc_prism.notna() & ch.depmap_id.notna() & ch.broad_id.notna()
        ].copy()
        if "passed_str_profiling" in ch:
            passed = ch.passed_str_profiling.astype("string").str.casefold().eq("true")
            ch = ch[passed.fillna(False)]
        chunks.append(ch.drop(columns=["auc"]))
    p = pd.concat(chunks, ignore_index=True)
    del chunks
    say(
        "DATA",
        f"GDSC usable={len(g):,} observations, {g.drug_id.nunique()} drug IDs, {g.model_id.nunique()} models; PRISM total={rows:,}, usable={len(p):,}, {p.broad_id.nunique()} compound IDs, {p.depmap_id.nunique()} cell lines",
    )
    # Avoid pretending multiple different assay/drug-batch entries are independent.
    drug_catalog = p.groupby("broad_id", as_index=False).agg(
        prism_name=("name", "first"), prism_n=("depmap_id", "nunique")
    )
    collisions = p.groupby("broad_id").name.nunique(dropna=True)
    if (collisions > 1).any():
        say(
            "CAUTION",
            f"{int((collisions > 1).sum())} PRISM broad IDs carry multiple names; quarantine those IDs",
        )
        drug_catalog = drug_catalog[
            ~drug_catalog.broad_id.isin(collisions[collisions > 1].index)
        ]
    gd = g.groupby("drug_id", as_index=False).agg(
        gdsc_name=("DRUG_NAME", "first"), gdsc_n=("model_id", "nunique")
    )
    gd["name_key"] = names(gd.gdsc_name)
    drug_catalog["name_key"] = names(drug_catalog.prism_name)
    candidates = gd.merge(drug_catalog, on="name_key", how="inner")
    candidates["match_status"] = "PROVISIONAL_NAME_ONLY"
    candidates.to_csv(OUT / "provisional_name_matches.csv", index=False)
    verified_file = OUT / "verified_compound_pairs.csv"
    if verified_file.exists():
        v = pd.read_csv(verified_file, dtype=str).fillna("")
        if not {"drug_id", "broad_id", "verification_evidence"}.issubset(v.columns):
            raise ValueError(
                "verified_compound_pairs.csv must contain drug_id,broad_id,verification_evidence"
            )
        v["drug_id"] = ident(v.drug_id)
        v["broad_id"] = ident(v.broad_id)
        v = v[
            (v.verification_evidence.str.strip() != "")
            & v.drug_id.isin(gd.drug_id)
            & v.broad_id.isin(drug_catalog.broad_id)
        ]
        v = v.drop_duplicates(["drug_id", "broad_id"])
    else:
        pd.DataFrame(columns=["drug_id", "broad_id", "verification_evidence"]).to_csv(
            verified_file, index=False
        )
        v = pd.DataFrame(columns=["drug_id", "broad_id", "verification_evidence"])
    # Distinct chemical identities in verified file must be checked by researcher.
    pairs = pd.concat(
        [
            candidates[
                ["drug_id", "broad_id", "gdsc_name", "prism_name", "match_status"]
            ],
            v.merge(gd[["drug_id", "gdsc_name"]], on="drug_id")
            .merge(drug_catalog[["broad_id", "prism_name"]], on="broad_id")
            .assign(match_status="USER_DOCUMENTED_VERIFICATION")[
                ["drug_id", "broad_id", "gdsc_name", "prism_name", "match_status"]
            ],
        ],
        ignore_index=True,
    )
    pairs = pairs.sort_values("match_status").drop_duplicates(
        ["drug_id", "broad_id"], keep="last"
    )
    pairs.to_csv(OUT / "candidate_pairs_for_analysis.csv", index=False)
    say(
        "MATCH",
        f"{len(candidates)} provisional name pairs; {len(v)} user-documented chemically verified pairs (if any).",
    )
    if not len(pairs):
        say(
            "CAUTION",
            "No candidate compounds matched; review exact names and chemical structures before evaluating.",
        )
    # Shared cell-line experimental concordance; compare ranks and signs, not absolute AUC values.
    result = []
    detail = []
    for pair in pairs.itertuples(index=False):
        a = g[g.drug_id == pair.drug_id][
            ["model_id", "depmap_id", "split", "auc_gdsc"]
        ].dropna(subset=["depmap_id"])
        b = p[p.broad_id == pair.broad_id][["depmap_id", "auc_prism"]].copy()
        if b.duplicated("depmap_id").any():
            # Multiple screens/replicates: no silent independent counting or cherry-picking.
            b = b.groupby("depmap_id", as_index=False).auc_prism.median()
        overlap = a.merge(b, on="depmap_id", validate="one_to_one")
        r = np.nan
        pv = np.nan
        if (
            len(overlap) >= 12
            and overlap.auc_gdsc.nunique() > 2
            and overlap.auc_prism.nunique() > 2
        ):
            test = spearmanr(overlap.auc_gdsc, overlap.auc_prism)
            r = float(test.statistic)
            pv = float(test.pvalue)
        result.append(
            dict(
                drug_id=pair.drug_id,
                broad_id=pair.broad_id,
                gdsc_name=pair.gdsc_name,
                prism_name=pair.prism_name,
                match_status=pair.match_status,
                n_gdsc=len(a),
                n_prism=b.depmap_id.nunique(),
                n_shared=len(overlap),
                n_gdsc_train=int(overlap.split.eq("train").sum()),
                n_gdsc_val=int(overlap.split.eq("val").sum()),
                n_gdsc_test=int(overlap.split.eq("test").sum()),
                spearman_auc=r,
                p=pv,
                platform_median_gdsc=float(overlap.auc_gdsc.median())
                if len(overlap)
                else np.nan,
                platform_median_prism=float(overlap.auc_prism.median())
                if len(overlap)
                else np.nan,
            )
        )
        if pair.drug_id in FOCUS:
            detail.extend(
                dict(
                    drug_id=pair.drug_id,
                    broad_id=pair.broad_id,
                    model_id=z.model_id,
                    depmap_id=z.depmap_id,
                    original_split=z.split,
                    gdsc_auc=z.auc_gdsc,
                    prism_auc=z.auc_prism,
                )
                for z in overlap.itertuples(index=False)
            )
    rr = pd.DataFrame(result)
    if not rr.empty:
        rr["q_fdr"] = np.nan
        tested = rr.p.notna()
        if tested.any():
            rr.loc[tested, "q_fdr"] = multipletests(
                rr.loc[tested, "p"].to_numpy(), method="fdr_bh"
            )[1]
        rr = rr.sort_values(["n_shared", "spearman_auc"], ascending=[False, False])
        rr.to_csv(OUT / "cross_platform_shared_line_concordance.csv", index=False)
    pd.DataFrame(
        detail,
        columns=[
            "drug_id",
            "broad_id",
            "model_id",
            "depmap_id",
            "original_split",
            "gdsc_auc",
            "prism_auc",
        ],
    ).to_csv(OUT / "focus_shared_cell_line_auc.csv", index=False)
    say(
        "COVERAGE",
        f"Candidate pairs with >=12 shared cell lines and rank correlations={int(rr.p.notna().sum()) if len(rr) else 0}",
    )
    if len(rr):
        say(
            "TOP CONCORDANCE",
            "\n"
            + pretty(
                rr[rr.n_shared >= 12].sort_values("spearman_auc", ascending=False),
                [
                    "gdsc_name",
                    "prism_name",
                    "match_status",
                    "n_shared",
                    "spearman_auc",
                    "q_fdr",
                ],
            ),
        )
        say(
            "BOTTOM CONCORDANCE",
            "\n"
            + pretty(
                rr[rr.n_shared >= 12].sort_values("spearman_auc"),
                [
                    "gdsc_name",
                    "prism_name",
                    "match_status",
                    "n_shared",
                    "spearman_auc",
                    "q_fdr",
                ],
            ),
        )
        say(
            "CASE STUDIES",
            "\n"
            + pretty(
                rr[rr.drug_id.isin(FOCUS)],
                [
                    "drug_id",
                    "gdsc_name",
                    "prism_name",
                    "match_status",
                    "n_shared",
                    "spearman_auc",
                    "q_fdr",
                ],
                25,
            ),
        )
    say(
        "MODEL TRANSFER",
        "Not performed: model transfer requires chemically verified pair mapping plus independently documented endpoint comparability and input feature coverage.",
    )
    say(
        "IMPORTANT",
        "Cross-platform correlations on shared cell lines are not independent-cell-line validation. Different PRISM/GDSC AUC protocols and curve fits may alter scale and ranks.",
    )
    say(
        "NEXT",
        "Inspect focus candidate pair chemistry and screen protocols; document matches in verified_compound_pairs.csv. Then plan leakage-free model-transfer evaluation and SLFN11/TP53 cross-platform hypotheses.",
    )
    say(
        "STATUS",
        f"SUCCESS; elapsed={(time.time() - START) / 60:.2f} min; results={OUT}",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as ex:
        say("ERROR", f"{type(ex).__name__}: {ex}")
        say("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 90 + "\n=== BEGIN PROJECT 5 PHASE 4A RESULTS ===")
        print("\n".join(REP))
        print("=== END PROJECT 5 PHASE 4A RESULTS ===\n" + "=" * 90)

# ----------------------------------------------------------------------------
# Corrected cross-platform biomarker analysis
# ----------------------------------------------------------------------------

from pathlib import Path
import time, traceback, warnings, re, zipfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
P2 = ROOT / "Results" / "Phase2"
P4A = ROOT / "Results" / "Phase4A"
OUT = ROOT / "Results" / "Phase4B"
OUT.mkdir(parents=True, exist_ok=True)
START = time.time()
LOG = []
FOCUS = {"1003": "Camptothecin", "1047": "Nutlin-3a (-)", "1060": "PD0325901"}
EXTRA = {
    "trametinib",
    "topotecan",
    "irinotecan",
    "selumetinib",
    "refametinib",
    "nutlin-3a",
    "nutlin-3a (-)",
}
MARKERS = ["SLFN11", "TP53", "BRAF", "KRAS", "NRAS", "IL6", "RSPO3", "ETV5", "CDKN2A"]


def say(tag, msg):
    LOG.append(f"[{tag}] {msg}")


def ids(v):
    return (
        v.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def normname(s):
    return (
        s.astype("string")
        .str.casefold()
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
    )


def f(prefix):
    files = list(DATA.glob(prefix + "*"))
    if len(files) != 1:
        raise FileNotFoundError(
            f"Expected one {prefix}* file in Data; found: {[a.name for a in files]}"
        )
    return files[0]


def table(x, n=20):
    return (
        x.head(n).to_string(index=False, float_format=lambda y: f"{y:.5g}")
        if len(x)
        else "None"
    )


def rho_p(x, y):
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 12 or np.unique(a[ok]).size < 3 or np.unique(b[ok]).size < 3:
        return int(ok.sum()), np.nan, np.nan
    res = spearmanr(a[ok], b[ok])
    return int(ok.sum()), float(res.statistic), float(res.pvalue)


def main():
    say(
        "DESIGN",
        "Paired GDSC2 vs PRISM associations on shared cell lines; no model retraining, no TEST outcome use and no direct AUC equivalence assumption.",
    )
    say(
        "DESIGN",
        "Compound-name matching is only provisional. This script DOES NOT automatically certify chemical structures or screen endpoint equivalence.",
    )
    candidates_path = P4A / "candidate_pairs_for_analysis.csv"
    if not candidates_path.exists():
        raise FileNotFoundError(
            "Run 05_cross_platform_validation.py from the beginning; Phase 4A candidate_pairs_for_analysis.csv is missing"
        )
    pairs = pd.read_csv(candidates_path, dtype={"drug_id": str, "broad_id": str})
    pairs["drug_id"] = ids(pairs.drug_id)
    pairs["broad_id"] = ids(pairs.broad_id)
    pairs["gdsc_key"] = normname(pairs.gdsc_name)
    chosen = pairs[pairs.drug_id.isin(FOCUS) | pairs.gdsc_key.isin(EXTRA)].copy()
    if chosen.empty:
        raise ValueError("No focus matched compound pairs in Phase4A candidate list")
    say(
        "MATCH",
        f"Phase4A candidate pairs={len(pairs)}; focus pairs={len(chosen)}; unique GDSC drug IDs={chosen.drug_id.nunique()}",
    )
    # Audit chemical identifiers from the PRISM side without claiming verification.
    ph = pd.read_csv(
        f("secondary-screen-dose-response-curve-parameters"), nrows=0
    ).columns.tolist()
    use = [
        x
        for x in [
            "broad_id",
            "depmap_id",
            "name",
            "auc",
            "smiles",
            "screen_id",
            "r2",
            "passed_str_profiling",
        ]
        if x in ph
    ]
    required = {"broad_id", "depmap_id", "name", "auc"}
    if not required.issubset(use):
        raise ValueError(f"Unexpected PRISM schema: {ph}")
    pchunks = []
    nrows = 0
    focus_b = set(chosen.broad_id)
    for part in pd.read_csv(
        f("secondary-screen-dose-response-curve-parameters"),
        usecols=use,
        chunksize=120000,
        low_memory=False,
    ):
        nrows += len(part)
        part["broad_id"] = ids(part.broad_id)
        part = part[part.broad_id.isin(focus_b)].copy()
        if "passed_str_profiling" in part.columns:
            part = part[
                part.passed_str_profiling.astype("string")
                .str.lower()
                .eq("true")
                .fillna(False)
            ]
        part["depmap_id"] = ids(part.depmap_id)
        part["prism_auc"] = pd.to_numeric(part.auc, errors="coerce")
        part = part[part.prism_auc.notna() & part.depmap_id.notna()].copy()
        if len(part):
            pchunks.append(part)
    if not pchunks:
        raise ValueError("No usable PRISM observations for focus drug candidates")
    prism = pd.concat(pchunks, ignore_index=True)
    del pchunks
    say(
        "PRISM",
        f"Rows scanned={nrows:,}; focus eligible rows={len(prism):,}; PRISM compounds={prism.broad_id.nunique()}",
    )
    summary = []
    ps = []
    for b, frame in prism.groupby("broad_id"):
        smiles = (
            sorted(set(frame.smiles.dropna().astype(str).str.strip()) - {"", "nan"})
            if "smiles" in frame
            else []
        )
        screens = int(frame.screen_id.nunique()) if "screen_id" in frame else np.nan
        summary.append(
            dict(
                broad_id=b,
                n_smiles=len(smiles),
                example_prism_smiles=smiles[0] if len(smiles) == 1 else "",
                n_screens=screens,
                n_cell_lines=frame.depmap_id.nunique(),
                n_rows=len(frame),
                n_multi_entry_cell_lines=int(frame.duplicated("depmap_id").sum()),
            )
        )
    chem = pd.DataFrame(summary)
    chem = chosen.merge(chem, on="broad_id", how="left")
    chem["verification_status"] = (
        "UNVERIFIED_NAME_MATCH_REQUIRES_GDSC_CHEMISTRY_AND_ENDPOINT_REVIEW"
    )
    chem["manual_evidence_required"] = (
        "Match stereo/tautomer/salt-aware chemical identifiers, dose ranges, protocol and active compound formulation"
    )
    chem.to_csv(OUT / "compound_identity_and_assay_audit.csv", index=False)
    say(
        "CHEMISTRY",
        f"PRISM identifiers present for {len(chem)} candidates; n_smiles=1 for {int((chem.n_smiles == 1).sum())}; NO chemical equivalence inferred from names or PRISM SMILES alone",
    )
    # Exact model IDs and pre-existing splits; never infer a sample association by model name.
    models = pd.read_csv(f("model_list_"), low_memory=False)
    if not {"model_id", "BROAD_ID"}.issubset(models):
        # Case-insensitive schema fallback, but never infer an unrelated ID field.
        bc = next(
            (c for c in models.columns if c.lower() in ("broad_id", "depmap_id")), None
        )
        if bc is None or "model_id" not in models:
            raise ValueError(f"Model identifiers missing: {list(models.columns)}")
    else:
        bc = "BROAD_ID"
    if "bc" not in locals():
        bc = "BROAD_ID"
    mapping = models[["model_id", bc]].copy().rename(columns={bc: "depmap_id"})
    mapping["model_id"] = ids(mapping.model_id)
    mapping["depmap_id"] = ids(mapping.depmap_id)
    mapping = mapping.dropna().drop_duplicates()
    # A model may have several recorded aliases. Never silently duplicate a model's
    # drug response across DepMap identifiers, or duplicate a DepMap result across
    # Sanger models. Exclude ambiguous mapping components for matched comparisons.
    n_raw = len(mapping)
    ambiguous_sanger = set(
        mapping.groupby("model_id").depmap_id.nunique().loc[lambda x: x > 1].index
    )
    ambiguous_depmap = set(
        mapping.groupby("depmap_id").model_id.nunique().loc[lambda x: x > 1].index
    )
    mapping = mapping[
        ~mapping.model_id.isin(ambiguous_sanger)
        & ~mapping.depmap_id.isin(ambiguous_depmap)
    ].copy()
    if mapping.duplicated("model_id").any() or mapping.duplicated("depmap_id").any():
        raise ValueError(
            "Mapping still non-bijective after excluding ambiguous aliases"
        )
    say(
        "ID AUDIT",
        f"Raw mapping pairs={n_raw:,}; ambiguous Sanger IDs={len(ambiguous_sanger)}; "
        f"ambiguous DepMap IDs={len(ambiguous_depmap)}; unique retained mappings={len(mapping):,}",
    )
    resp = pd.read_excel(
        f("GDSC2_fitted_dose_response_"), usecols=["SANGER_MODEL_ID", "DRUG_ID", "AUC"]
    )
    resp["model_id"] = ids(resp.SANGER_MODEL_ID)
    resp["drug_id"] = ids(resp.DRUG_ID)
    resp["gdsc_auc"] = pd.to_numeric(resp.AUC, errors="coerce")
    resp = resp[resp.drug_id.isin(chosen.drug_id) & resp.gdsc_auc.between(0, 1)][
        ["model_id", "drug_id", "gdsc_auc"]
    ]
    if resp.duplicated(["model_id", "drug_id"]).any():
        raise ValueError("Duplicate GDSC drug×model entries")
    splits = pd.read_csv(P2 / "model_splits.csv")
    splits["model_id"] = ids(splits.model_id)
    if splits.model_id.duplicated().any():
        raise ValueError("Duplicate model IDs in splits")
    # Selected RNA matrix may not contain all relevant genes; explicitly report unavailable markers.
    rna = pd.read_csv(
        P2 / "expression_selected_genes.csv.gz", compression="gzip", low_memory=False
    )
    rna["model_id"] = ids(rna.model_id)
    present = [x for x in MARKERS if x in rna.columns]
    absent = [x for x in MARKERS if x not in rna.columns]
    say(
        "RNA",
        f"Phase2 selected RNA matrix contains markers={present}; absent={absent} (absent is NOT expression-negative)",
    )
    rna = rna[["model_id"] + present].drop_duplicates("model_id")
    # Add only *recorded* mutation events, not inferred WT. Summary is not whole exome.
    with zipfile.ZipFile(f("mutations_summary_")) as z:
        poss = [
            x
            for x in z.namelist()
            if x.lower().endswith(".csv") and "mutations_summary_" in x.lower()
        ]
        if len(poss) != 1:
            raise ValueError(f"Cannot uniquely identify mutation summary: {poss}")
        with z.open(poss[0]) as q:
            mut = pd.read_csv(q, usecols=["model_id", "gene_symbol"], low_memory=False)
    mut["model_id"] = ids(mut.model_id)
    mut["gene_symbol"] = mut.gene_symbol.astype("string").str.strip()
    relevant = {"TP53", "BRAF", "KRAS", "NRAS"}
    mut = mut[mut.gene_symbol.isin(relevant)].drop_duplicates(
        ["model_id", "gene_symbol"]
    )
    mutset = {
        gene: set(mut.loc[mut.gene_symbol == gene, "model_id"]) for gene in relevant
    }
    # Critical: restrict mutation comparisons to metadata-confirmed WES-assayed cell lines.
    av = pd.read_csv(f("model_dataset_availability_"), low_memory=False)
    av["model_id"] = ids(av.model_id)
    coverage_cols = [
        x
        for x in av.columns
        if x in ("Mutation Sanger Cell Lines (WES)", "Mutation Broad Cell Lines (WES)")
    ]
    if not coverage_cols:
        raise ValueError("Missing WES assay coverage metadata")
    covered = pd.Series(False, index=av.index)
    for c in coverage_cols:
        covered |= av[c].astype("string").str.lower().eq("true").fillna(False)
    covered_ids = set(av.loc[covered, "model_id"])
    say(
        "MUTATION",
        f"WES coverage-flagged models={len(covered_ids)}; summary records are event-present vs no recorded event, NOT confirmed mutant vs WT",
    )
    d = pd.DataFrame()
    details = []
    results = []
    comparisons = []
    for row in chosen.itertuples(index=False):
        sub = prism[prism.broad_id == row.broad_id].copy()
        # preserve information about repeated screens, prevent pseudoreplication.
        repeated = int(sub.duplicated("depmap_id").sum())
        sub = sub.groupby("depmap_id", as_index=False).prism_auc.median()
        gg = resp[resp.drug_id == row.drug_id].merge(
            mapping, on="model_id", how="inner", validate="many_to_one"
        )
        gg = gg.merge(
            splits[["model_id", "split"]],
            on="model_id",
            how="left",
            validate="one_to_one",
        )
        if gg.duplicated("depmap_id").any():
            raise ValueError(
                f"Unexpected duplicate mapped DepMap ID for GDSC drug {row.drug_id}"
            )
        joined = gg.merge(sub, on="depmap_id", validate="one_to_one")
        joined = joined.merge(rna, on="model_id", how="left", validate="one_to_one")
        if len(joined):
            joined["drug_id"] = row.drug_id
            joined["broad_id"] = row.broad_id
            details.append(joined)
        n, rs, pp = rho_p(joined.gdsc_auc, joined.prism_auc)
        comparisons.append(
            dict(
                drug_id=row.drug_id,
                broad_id=row.broad_id,
                drug=row.gdsc_name,
                n_shared=n,
                spearman_rho=rs,
                p=pp,
                duplicated_prism_rows=repeated,
                original_train=int(joined.split.eq("train").sum()),
                original_val=int(joined.split.eq("val").sum()),
                original_test=int(joined.split.eq("test").sum()),
            )
        )
        for marker in present:
            y = pd.to_numeric(joined[marker], errors="coerce")
            for platform, endpoint in [("GDSC2", "gdsc_auc"), ("PRISM", "prism_auc")]:
                nn, rr, px = rho_p(y, joined[endpoint])
                results.append(
                    dict(
                        drug_id=row.drug_id,
                        broad_id=row.broad_id,
                        drug=row.gdsc_name,
                        marker=marker,
                        modality="RNA",
                        platform=platform,
                        n=nn,
                        effect=rr,
                        p=px,
                        analysis="SPEARMAN",
                    )
                )
        if len(joined):
            eligible = joined[joined.model_id.isin(covered_ids)]
            for gene in relevant:
                flag = eligible.model_id.isin(mutset[gene])
                for platform, endpoint in [
                    ("GDSC2", "gdsc_auc"),
                    ("PRISM", "prism_auc"),
                ]:
                    yy = eligible[endpoint]
                    a = yy[flag].dropna()
                    b = yy[~flag].dropna()
                    pv = (
                        float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
                        if len(a) >= 5 and len(b) >= 5
                        else np.nan
                    )
                    results.append(
                        dict(
                            drug_id=row.drug_id,
                            broad_id=row.broad_id,
                            drug=row.gdsc_name,
                            marker=gene,
                            modality="RECORDED_MUTATION_EVENT",
                            platform=platform,
                            n=len(a) + len(b),
                            n_event=len(a),
                            n_no_record=len(b),
                            effect=float(a.median() - b.median())
                            if len(a) and len(b)
                            else np.nan,
                            p=pv,
                            analysis="MWU_MEDIAN_DIFFERENCE",
                        )
                    )
    concord = pd.DataFrame(comparisons)
    concord.to_csv(OUT / "focus_crossplatform_concordance.csv", index=False)
    associations = pd.DataFrame(results)
    if len(associations):
        associations["q_fdr"] = np.nan
        ok = associations.p.notna()
        if ok.any():
            associations.loc[ok, "q_fdr"] = multipletests(
                associations.loc[ok, "p"], method="fdr_bh"
            )[1]
        associations["direction"] = np.where(
            associations.effect > 0,
            "higher_AUC",
            np.where(associations.effect < 0, "lower_AUC", "undefined"),
        )
        associations.to_csv(
            OUT / "focus_crossplatform_biomarker_associations.csv", index=False
        )
    if details:
        pd.concat(details, ignore_index=True).to_csv(
            OUT / "focus_shared_cellline_molecular_table.csv", index=False
        )
    say(
        "CONCORDANCE",
        "\n" + table(concord.sort_values("spearman_rho", ascending=False), 25),
    )
    say(
        "BIOMARKERS",
        f"Number of association tests={len(associations)}, usable={int(associations.p.notna().sum())}, exploratory q<0.05={int((associations.q_fdr < 0.05).sum())}",
    )
    top = associations[associations.p.notna()].sort_values("q_fdr")
    say(
        "TOP ASSOCIATIONS",
        "\n"
        + table(
            top[
                ["drug", "marker", "modality", "platform", "n", "effect", "p", "q_fdr"]
            ],
            25,
        ),
    )
    replicated = []
    for key, group in associations.groupby(
        ["drug_id", "broad_id", "marker", "modality", "analysis"]
    ):
        if set(group.platform) == {"PRISM", "GDSC2"}:
            a = group.set_index("platform")
            x = a.loc["GDSC2"]
            y = a.loc["PRISM"]
            replicated.append(
                dict(
                    drug_id=key[0],
                    broad_id=key[1],
                    marker=key[2],
                    modality=key[3],
                    n_gdsc=int(x["n"]),
                    n_prism=int(y["n"]),
                    effect_gdsc=x["effect"],
                    effect_prism=y["effect"],
                    q_gdsc=x.q_fdr,
                    q_prism=y.q_fdr,
                    same_direction=bool(
                        np.isfinite(x["effect"])
                        and np.isfinite(y["effect"])
                        and np.sign(x["effect"]) == np.sign(y["effect"])
                    ),
                    both_q_lt_05=bool(x.q_fdr < 0.05 and y.q_fdr < 0.05),
                )
            )
    rep = pd.DataFrame(replicated)
    rep.to_csv(OUT / "paired_biomarker_direction_consistency.csv", index=False)
    if len(rep):
        say(
            "REPLICATION",
            f"Paired associations={len(rep)}; matching directional signs={rep.same_direction.sum()}; both exploratory FDR<.05={rep.both_q_lt_05.sum()}",
        )
        say("BOTH SIGNIFICANT", "\n" + table(rep[rep.both_q_lt_05], 25))
    # Identify the precise verification gate for subsequent frozen-model deployment.
    gate = chem[
        [
            "drug_id",
            "broad_id",
            "gdsc_name",
            "prism_name",
            "verification_status",
            "n_smiles",
            "n_screens",
            "n_multi_entry_cell_lines",
        ]
    ].copy()
    gate = gate.merge(
        concord[
            [
                "drug_id",
                "broad_id",
                "n_shared",
                "original_train",
                "original_val",
                "original_test",
            ]
        ],
        on=["drug_id", "broad_id"],
        how="left",
    )
    gate["model_transfer_eligible_now"] = False
    gate["blocker"] = (
        "Chemical identity unverified; PRISM and GDSC AUC endpoint/dose comparability not established; matched lines are NOT unseen models"
    )
    gate.to_csv(OUT / "model_transfer_evidence_gate.csv", index=False)
    say(
        "GATE",
        f"Model transfer: 0 certified-eligible pairs; saved chemical/assay audit and evidence-gate CSV for review.",
    )
    say(
        "CAUTION",
        "Both-platform biomarker associations share the SAME cell lines; directional agreement is cross-assay replication, not independent biological replication.",
    )
    say(
        "CAUTION",
        "PRISM multiple response entries per line are summarised by median as a pragmatic exploratory step; inspect screen_id and replicate definitions before confirmatory analysis.",
    )
    say(
        "CAUTION",
        "RNA from preselected 256-gene model matrix only; a key biomarker absent from it is UNTESTED, not nonsignificant. RNA expression may reflect Broad/Sanger technical source.",
    )
    say(
        "NEXT",
        "Inspect identity audit; check GDSC structures/identifiers and PRISM compound SMILES, salt/stereochemistry and screen ranges. Use compatible frozen-model evaluation only if verification gates pass.",
    )
    say(
        "STATUS", f"SUCCESS; elapsed={(time.time() - START) / 60:.2f} min; output={OUT}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as ex:
        say("ERROR", f"{type(ex).__name__}: {ex}")
        say("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 95 + "\n=== BEGIN PROJECT 5 PHASE 4B RESULTS ===")
        print("\n".join(LOG))
        print("=== END PROJECT 5 PHASE 4B RESULTS ===\n" + "=" * 95)

# ----------------------------------------------------------------------------
# Adjusted association robustness
# ----------------------------------------------------------------------------

from pathlib import Path
import zipfile, time, traceback, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, t as student_t
from statsmodels.stats.multitest import multipletests

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
P4B = ROOT / "Results" / "Phase4B"
OUT = ROOT / "Results" / "Phase4C"
OUT.mkdir(parents=True, exist_ok=True)
LOG = []
START = time.time()
GENES = ["SLFN11", "IL6", "RSPO3", "ETV5", "CDKN2A", "TP53", "BRAF", "KRAS", "NRAS"]
FOCUS_NAMES = {
    "camptothecin",
    "irinotecan",
    "topotecan",
    "trametinib",
    "selumetinib",
    "refametinib",
}
MIN_GROUP = 12
MIN_OVERALL = 35


def say(k, x):
    LOG.append(f"[{k}] {x}")


def ident(s):
    return (
        s.astype("string").str.strip().str.upper().str.replace(r"\.0$", "", regex=True)
    )


def file(prefix):
    a = list(DATA.glob(prefix + "*"))
    if len(a) != 1:
        raise FileNotFoundError(
            f"Expected one Data/{prefix}* file, found {[p.name for p in a]}"
        )
    return a[0]


def show(df, n=16):
    return (
        "No eligible data"
        if df.empty
        else df.head(n).to_string(index=False, float_format=lambda v: f"{v:.4g}")
    )


def safe_spearman(x, y, min_n=12):
    a = np.asarray(x, float)
    b = np.asarray(y, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < min_n or len(np.unique(a[ok])) < 3 or len(np.unique(b[ok])) < 3:
        return int(ok.sum()), np.nan, np.nan
    z = spearmanr(a[ok], b[ok])
    return int(ok.sum()), float(z.statistic), float(z.pvalue)


def adjusted_rank_assoc(df, gene, outcome, adjust):
    """Rank residual correlation adjusted for categorical cancer and RNA source.
    p-value uses approximate partial-correlation t reference; exploratory only.
    """
    d = df[[gene, outcome] + adjust].dropna().copy()
    n = len(d)
    if n < MIN_OVERALL or d[gene].nunique() < 4 or d[outcome].nunique() < 4:
        return n, np.nan, np.nan, 0
    x = rankdata(pd.to_numeric(d[gene]).to_numpy(float))
    y = rankdata(pd.to_numeric(d[outcome]).to_numpy(float))
    if adjust:
        dummy = pd.get_dummies(d[adjust].astype(str), drop_first=True, dtype=float)
        Z = np.column_stack([np.ones(n), dummy.to_numpy(float)])
    else:
        Z = np.ones((n, 1))
    effective = np.linalg.matrix_rank(Z)
    if n - effective - 1 < 6:
        return n, np.nan, np.nan, effective
    xx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    yy = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    den = np.linalg.norm(xx) * np.linalg.norm(yy)
    if den < 1e-12:
        return n, np.nan, np.nan, effective
    r = float(np.clip(xx @ yy / den, -1, 1))
    dfree = n - effective - 1
    t = abs(r) * np.sqrt(dfree / max(1 - r * r, 1e-12))
    return n, r, float(2 * student_t.sf(t, dfree)), effective


def load_full_rna(models):
    """Extract ONLY requested genes and cell lines from full RNA ZIP, with source metadata.
    Three metadata rows are model_name, data_source, gene_symbol (confirmed in Phase2D).
    """
    with zipfile.ZipFile(file("rnaseq_merged_")) as z:
        cand = [
            n
            for n in z.namelist()
            if "rnaseq_merged_rsem_tpm_" in Path(n).name and n.lower().endswith(".csv")
        ]
        if len(cand) != 1:
            raise ValueError(f"Unexpected RNA archive members: {cand}")
        with z.open(cand[0]) as h:
            head = pd.read_csv(h, nrows=0).columns.tolist()
        present = sorted(set(head) & set(models))
        use = head[:3] + present
        if not present:
            raise ValueError("No matching Sanger model IDs in RNA source")
        # Restrict columns at parse time. Retaining all genes avoids unsupported gene absence.
        with z.open(cand[0]) as h:
            r = pd.read_csv(h, usecols=use, dtype="string", low_memory=False)
    if list(r.iloc[:3, 0]) != ["model_name", "data_source", "gene_symbol"]:
        raise ValueError(
            f"RNA metadata rows unexpectedly differ: {r.iloc[:3, 0].tolist()}"
        )
    sources = r.loc[1, present].astype(str).to_dict()
    # IMPORTANT: the first column contains gene SYMBOLS in data rows.
    # 'Unnamed: 2' contains numeric gene IDs; selecting symbols there silently gives zero matches.
    gene_symbols = r.iloc[3:, 0].astype("string").str.strip().to_numpy()
    sel = np.isin(gene_symbols, GENES)
    part = r.iloc[3:].loc[sel, ["model_id"] + present].copy()
    part = part.rename(columns={"model_id": "gene"})
    if part.gene.duplicated().any():
        dup = part.loc[part.gene.duplicated(False), "gene"].tolist()
        say(
            "RNA",
            f"Duplicate requested symbols excluded (no arbitrary aggregation): {dup}",
        )
        part = part[~part.gene.duplicated(False)]
    if part.empty:
        raise ValueError(
            "No prespecified gene symbols matched full RNA. Check the model_id gene-symbol column and ZIP schema."
        )
    long = part.melt(id_vars="gene", var_name="model_id", value_name="tpm")
    long["tpm"] = pd.to_numeric(long.tpm, errors="coerce")
    long.loc[long.tpm < 0, "tpm"] = np.nan
    long["expression"] = np.log2(long.tpm + 1)
    wide = long.pivot(
        index="model_id", columns="gene", values="expression"
    ).reset_index()
    wide["rna_source"] = wide.model_id.map(sources)
    say(
        "RNA",
        f"Full RNA source matched models={len(present)}; requested genes present={sorted(part.gene.dropna().unique())}; sources={wide.rna_source.value_counts(dropna=False).to_dict()}",
    )
    return wide


def main():
    say(
        "DESIGN",
        "Exploratory GDSC2-PRISM biomarker association robustness; no new AI training, model transfer or TEST evaluation.",
    )
    say(
        "DESIGN",
        "Uses full RNA for prespecified markers, partial rank associations adjusted for cancer type and RNA source.",
    )
    say(
        "DESIGN",
        "GDSC2/PRISM matched cell lines are the SAME biological models; cross-assay agreement is not unseen-cell-line validation.",
    )
    say(
        "DESIGN",
        "All chemical name matches remain PROVISIONAL. PRISM AUC is not numerically equivalent to GDSC2 AUC.",
    )
    shared = P4B / "focus_shared_cellline_molecular_table.csv"
    aud = P4B / "compound_identity_and_assay_audit.csv"
    if not shared.exists() or not aud.exists():
        raise FileNotFoundError("Phase4B files missing. Run fixed Phase4B first.")
    a = pd.read_csv(aud, dtype={"drug_id": str, "broad_id": str})
    say(
        "IDENTITY",
        f"Audited candidate pairs={len(a)}; all status={a.verification_status.value_counts(dropna=False).to_dict()}",
    )
    d = pd.read_csv(
        shared,
        dtype={"drug_id": str, "broad_id": str, "model_id": str, "depmap_id": str},
        low_memory=False,
    )
    for c in ["model_id", "depmap_id", "drug_id", "broad_id"]:
        d[c] = ident(d[c])
    if d.duplicated(["drug_id", "broad_id", "model_id"]).any():
        raise ValueError("Unexpected repeated pair × model rows in Phase4B output")
    meta = pd.read_csv(
        file("model_list_"),
        usecols=["model_id", "cancer_type", "tissue_status", "model_type"],
        low_memory=False,
    )
    meta.model_id = ident(meta.model_id)
    if meta.model_id.duplicated().any():
        raise ValueError("Duplicate model IDs in Cell Model Passports metadata")
    d = d.merge(meta, on="model_id", how="left", validate="many_to_one")
    d = d[
        (d.tissue_status.isin(["Tumour", "Metastasis"]))
        & d.model_type.eq("Cell Line")
        & d.cancer_type.notna()
    ].copy()
    if d.empty:
        raise ValueError("No tumour cell lines after metadata filtering")
    r = load_full_rna(d.model_id.unique())
    d = d.drop(columns=[c for c in GENES if c in d.columns], errors="ignore").merge(
        r, on="model_id", how="left", validate="many_to_one"
    )
    pairs = a[["drug_id", "broad_id", "gdsc_name", "prism_name"]].drop_duplicates(
        ["drug_id", "broad_id"]
    )
    d = d.merge(pairs, on=["drug_id", "broad_id"], how="left", validate="many_to_one")
    d.to_csv(OUT / "matched_cells_full_marker_matrix.csv", index=False)
    say(
        "COHORT",
        f"Matched tumour pair-model records={len(d)}; unique cell lines={d.model_id.nunique()}; pairs={d.groupby(['drug_id', 'broad_id']).ngroups}",
    )
    rows = []
    within = []
    for (gid, bid), f in d.groupby(["drug_id", "broad_id"], sort=False):
        name = str(f.gdsc_name.iloc[0])
        nn = f.model_id.nunique()
        for gene in GENES:
            if gene not in f.columns:
                continue
            for platform, y in [("GDSC2", "gdsc_auc"), ("PRISM", "prism_auc")]:
                if y not in f:
                    raise KeyError(f"Missing required Phase4B outcome {y}")
                n, rho, p = safe_spearman(f[gene], f[y])
                nz, partial, p_adj, rank = adjusted_rank_assoc(
                    f, gene, y, ["cancer_type", "rna_source"]
                )
                rows.append(
                    dict(
                        drug_id=gid,
                        broad_id=bid,
                        drug=name,
                        gene=gene,
                        platform=platform,
                        n=n,
                        unadjusted_rho=rho,
                        unadjusted_p=p,
                        adjusted_n=nz,
                        adjusted_rho=partial,
                        adjusted_p=p_adj,
                        covariate_rank=rank,
                        n_cancer_types=f.cancer_type.nunique(),
                    )
                )
                if gene == "SLFN11" and name.lower() in FOCUS_NAMES:
                    for cancer, g in f.groupby("cancer_type"):
                        if len(g) < MIN_GROUP:
                            continue
                        gn, gr, gp = safe_spearman(g[gene], g[y], MIN_GROUP)
                        within.append(
                            dict(
                                drug_id=gid,
                                broad_id=bid,
                                drug=name,
                                gene=gene,
                                platform=platform,
                                cancer_type=cancer,
                                n=gn,
                                rho=gr,
                                p=gp,
                            )
                        )
    res = pd.DataFrame(
        rows,
        columns=[
            "drug_id",
            "broad_id",
            "drug",
            "gene",
            "platform",
            "n",
            "unadjusted_rho",
            "unadjusted_p",
            "adjusted_n",
            "adjusted_rho",
            "adjusted_p",
            "covariate_rank",
            "n_cancer_types",
        ],
    )
    res["q_adjusted"] = np.nan
    mask = res.adjusted_p.notna()
    if mask.any():
        res.loc[mask, "q_adjusted"] = multipletests(
            res.loc[mask, "adjusted_p"], method="fdr_bh"
        )[1]
    res.to_csv(OUT / "cancer_source_adjusted_associations.csv", index=False)
    within = pd.DataFrame(within)
    if len(within):
        within["q_fdr"] = np.nan
        ok = within.p.notna()
        if ok.any():
            within.loc[ok, "q_fdr"] = multipletests(
                within.loc[ok, "p"], method="fdr_bh"
            )[1]
        within.to_csv(OUT / "within_cancer_SLFN11_associations.csv", index=False)
    say(
        "ADJUSTMENT",
        f"Assay-drug-gene contrasts={len(res)}; testable after adjustment={mask.sum()}; exploratory adjusted FDR<0.05={(res.q_adjusted < 0.05).sum()}",
    )
    target = res[res.gene.eq("SLFN11")].sort_values(["drug", "broad_id", "platform"])
    say(
        "SLFN11 ALL",
        "\n"
        + show(
            target[
                [
                    "drug",
                    "drug_id",
                    "broad_id",
                    "platform",
                    "n",
                    "unadjusted_rho",
                    "adjusted_rho",
                    "q_adjusted",
                ]
            ],
            26,
        ),
    )
    interesting = res[res.gene.isin(["IL6", "RSPO3", "ETV5"])].sort_values(
        "q_adjusted", na_position="last"
    )
    say(
        "IL6/RSPO3/ETV5",
        "\n"
        + show(
            interesting[
                ["drug", "gene", "platform", "n", "adjusted_rho", "q_adjusted"]
            ],
            18,
        ),
    )
    sig = res[(res.q_adjusted < 0.05)].sort_values("q_adjusted")
    say(
        "OTHER SUPPORTED ASSOCIATIONS",
        "\n"
        + show(
            sig[["drug", "gene", "platform", "n", "adjusted_rho", "q_adjusted"]], 20
        ),
    )
    if len(within):
        say(
            "WITHIN CANCER",
            f"Tested={within.p.notna().sum()}, FDR<.05={(within.q_fdr < 0.05).sum()}; significant examples:\n"
            + show(
                within[within.q_fdr < 0.05].sort_values("q_fdr")[
                    ["drug", "platform", "cancer_type", "n", "rho", "q_fdr"]
                ],
                18,
            ),
        )
    # Pair platform results; do not count distinct BRD assay formulations as independent biological replication.
    both = res.pivot_table(
        index=["drug_id", "broad_id", "drug", "gene"],
        columns="platform",
        values=["adjusted_rho", "q_adjusted"],
        aggfunc="first",
    )
    both.columns = ["_".join(map(str, c)) for c in both.columns]
    both = both.reset_index()
    need = [
        "adjusted_rho_GDSC2",
        "adjusted_rho_PRISM",
        "q_adjusted_GDSC2",
        "q_adjusted_PRISM",
    ]
    if all(c in both for c in need):
        both["same_direction"] = np.sign(both.adjusted_rho_GDSC2) == np.sign(
            both.adjusted_rho_PRISM
        )
        both["both_fdr_lt05"] = (both.q_adjusted_GDSC2 < 0.05) & (
            both.q_adjusted_PRISM < 0.05
        )
        both["association_support"] = "UNADJUSTED_OR_PLATFORM_SPECIFIC"
        both.loc[both.same_direction & both.both_fdr_lt05, "association_support"] = (
            "BOTH_ASSAYS_EXPLORATORY_FDR"
        )
        both.to_csv(OUT / "paired_adjusted_association_evidence.csv", index=False)
        say(
            "PAIRED",
            f"Pair-marker entries={len(both)}; matching adjusted directions={both.same_direction.sum()}; both-adjusted-FDR<.05={both.both_fdr_lt05.sum()}",
        )
        say(
            "PAIRED TOP",
            "\n"
            + show(
                both[both.same_direction & both.both_fdr_lt05][
                    [
                        "drug",
                        "gene",
                        "adjusted_rho_GDSC2",
                        "adjusted_rho_PRISM",
                        "q_adjusted_GDSC2",
                        "q_adjusted_PRISM",
                    ]
                ],
                20,
            ),
        )
    # Repurposing gate: indications are not interchangeable with assay sensitivity.
    gate = pairs.copy()
    gate["chemical_identity_verified"] = False
    gate["assay_equivalence_verified"] = False
    gate["indication_novelty_checked"] = False
    gate["cross_platform_unseen_cellline_prediction"] = False
    gate["status"] = "NOT_READY_FOR_REPURPOSING_CLAIM"
    gate["next_requirements"] = (
        "Indication/approval review; compound structure and dose/protocol audit; external clinical/translational context; genuine unseen-model tests"
    )
    gate.to_csv(OUT / "repurposing_evidence_gates.csv", index=False)
    say(
        "REPURPOSING",
        f"Candidates from checked drug pairs={len(gate)}; evidence-qualified repurposing claims=0 (by design: no indication/chemistry gate passed).",
    )
    say(
        "LIMITS",
        "Approximate partial-Spearman p-values, source/cancer confounding, shared cell lines, correlated drug pairs and reused exploratory cohorts. Source adjustment cannot fully remove batch confounding.",
    )
    say(
        "LIMITS",
        "Restricted to name-provisional selected compounds. No structural or dose-range verification, no external AI model-transfer, and no clinical conclusions.",
    )
    say(
        "NEXT",
        "Review cancer/source-adjusted SLFN11 results and per-cancer n; verify chemical identities and assay protocols before frozen-model transfer; independently review indications before naming repurposing candidates.",
    )
    say(
        "STATUS",
        f"SUCCESS; elapsed={(time.time() - START) / 60:.2f} min; outputs={OUT}",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        say("ERROR", f"{type(exc).__name__}: {exc}")
        say("TRACEBACK", traceback.format_exc())
        raise
    finally:
        print("\n" + "=" * 100 + "\n=== BEGIN PROJECT 5 PHASE 4C RESULTS ===")
        print("\n".join(LOG))
        print("=== END PROJECT 5 PHASE 4C RESULTS ===\n" + "=" * 100)
