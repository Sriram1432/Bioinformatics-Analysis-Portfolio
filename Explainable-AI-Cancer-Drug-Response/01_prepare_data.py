"""Validate and harmonise the pharmacogenomic source data.

Audits GDSC2 responses, model identifiers and multi-omics availability without treating cell lines as patients.

Run from the repository root with ``python 01_prepare_data.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# Source inventory and response-data integrity
# ----------------------------------------------------------------------------

import re, json, math, sys, traceback, warnings, zipfile, io
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase1"
OUT.mkdir(parents=True, exist_ok=True)
REPORT = []
ISSUES = []
CATALOG = []
TABLES = {}
MAX_INSPECTION_ROWS = 1200  # sample only: large omics files
PRISM_CHUNK = 150_000


def note(section, key, value):
    s = str(value).replace("\n", " | ")
    REPORT.append((section, key, s[:1700]))


def warn(stage, error):
    msg = f"{stage}: {str(error)[:450]}"
    ISSUES.append(msg)
    note("WARNINGS", stage, msg)


def tidy(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).lower()).strip("_")


def clean_id(s):
    x = s.astype("string").str.strip().str.upper()
    x = x.str.replace(r"\.0$", "", regex=True)
    return x.mask(x.isin(["", "NAN", "NONE", "<NA>", "NULL"]))


def match_col(df, exact=(), words=(), excludes=()):
    d = {tidy(c): c for c in df.columns}
    for c in exact:
        if tidy(c) in d:
            return d[tidy(c)]
    if words:
        for c in df.columns:
            t = tidy(c)
            if all(w in t for w in words) and not any(z in t for z in excludes):
                return c
    return None


def safe_top(s, n=12):
    return "; ".join(
        f"{str(k)[:38]}={int(v):,}"
        for k, v in s.astype("string")
        .fillna("Missing")
        .value_counts(dropna=False)
        .head(n)
        .items()
    )


def qstats(s):
    x = pd.to_numeric(s, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if not len(x):
        return "no numeric measurements"
    q = x.quantile([0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1])
    return f"n={len(x):,}; median={q.loc[0.5]:.4g}; q05={q.loc[0.05]:.4g}; q95={q.loc[0.95]:.4g}; min={q.loc[0]:.4g}; max={q.loc[1]:.4g}"


def csv_read(path, nrows=None, chunksize=None):
    kw = dict(low_memory=False)
    if nrows is not None:
        kw["nrows"] = nrows
    if chunksize is not None:
        kw["chunksize"] = chunksize
    try:
        return pd.read_csv(path, **kw)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="latin1", **kw)


def peek(path, n=MAX_INSPECTION_ROWS):
    suf = path.suffix.lower()
    if suf in (".xlsx", ".xls"):
        return pd.read_excel(
            path, nrows=n, engine="openpyxl" if suf == ".xlsx" else None
        )
    if suf in (".csv", ".tsv", ".txt"):
        try:
            return pd.read_csv(
                path, nrows=n, sep="\t" if suf == ".tsv" else None, engine="python"
            )
        except UnicodeDecodeError:
            return pd.read_csv(
                path,
                nrows=n,
                sep="\t" if suf == ".tsv" else None,
                engine="python",
                encoding="latin1",
            )
    if suf == ".zip":
        with zipfile.ZipFile(path) as z:
            items = [m for m in z.infolist() if not m.is_dir()]
            choices = [
                m
                for m in items
                if Path(m.filename).suffix.lower() in (".csv", ".tsv", ".txt", ".xlsx")
            ]
            if not choices:
                return None
            m = sorted(
                choices, key=lambda x: (x.file_size > 2_000_000_000, -x.file_size)
            )[0]
            with z.open(m) as fh:
                if m.filename.lower().endswith(".xlsx"):
                    return pd.read_excel(io.BytesIO(fh.read()), nrows=n)
                buf = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
                return pd.read_csv(
                    buf,
                    nrows=n,
                    sep="\t" if m.filename.endswith(".tsv") else None,
                    engine="python",
                )
    return None


def save(df, name):
    df.to_csv(OUT / name, index=False)


def inventory():
    note("SETUP", "Project folder", str(ROOT))
    note("SETUP", "Python / Pandas", f"{sys.version.split()[0]} / {pd.__version__}")
    if not DATA.exists():
        warn("DATA DIRECTORY", f"Folder missing: {DATA}")
        return
    for p in sorted(DATA.iterdir()):
        if not p.is_file():
            continue
        row = dict(
            file=p.name,
            size_mb=round(p.stat().st_size / 1024**2, 2),
            members="",
            columns="",
            inspected_rows=0,
        )
        try:
            if p.suffix.lower() == ".zip":
                with zipfile.ZipFile(p) as z:
                    members = [m for m in z.infolist() if not m.is_dir()]
                    row["members"] = len(members)
                    member_names = [
                        f"{m.filename} ({m.file_size / 1024**2:.1f}MB)"
                        for m in sorted(members, key=lambda x: -x.file_size)[:12]
                    ]
                    note("ZIP CONTENTS", p.name, "; ".join(member_names))
            if p.suffix.lower() in (".csv", ".tsv", ".xlsx", ".xls", ".zip", ".txt"):
                a = peek(p, 40)
                if a is not None:
                    row["columns"] = " | ".join(map(str, a.columns[:45]))
                    row["inspected_rows"] = len(a)
                    note("SOURCE COLUMNS", p.name, row["columns"])
        except Exception as exc:
            warn("inventory " + p.name, exc)
        CATALOG.append(row)
    if CATALOG:
        save(pd.DataFrame(CATALOG), "file_inventory.csv")
    note("SETUP", "Source files", len(CATALOG))


def getfile(fragment):
    files = [
        p for p in DATA.iterdir() if p.is_file() and fragment.lower() in p.name.lower()
    ]
    return max(files, key=lambda p: p.stat().st_mtime) if files else None


def read_small(fragment, required=False):
    p = getfile(fragment)
    if not p:
        if required:
            warn("missing " + fragment, "not found")
        return None
    try:
        if p.suffix.lower() == ".zip":
            with zipfile.ZipFile(p) as z:
                cands = [n for n in z.namelist() if n.lower().endswith(".csv")]
                cands = sorted(cands, key=lambda name: z.getinfo(name).file_size)
                if not cands:
                    return peek(p)
                with z.open(cands[0]) as f:
                    return pd.read_csv(f, low_memory=False)
        if p.suffix.lower() == ".csv":
            return csv_read(p)
        if p.suffix.lower() == ".xlsx":
            return pd.read_excel(p)
        return peek(p)
    except Exception as exc:
        warn("reading " + fragment, exc)
        return None


def model_meta():
    df = read_small("model_list_", True)
    if df is None:
        return {}
    note("MODELS", "Model rows", len(df))
    note("MODELS", "Columns", " | ".join(map(str, df.columns)))
    for c in df.columns:
        t = tidy(c)
        if any(
            w in t
            for w in (
                "tissue",
                "cancer_type",
                "model_type",
                "msi_status",
                "growth",
                "cosmic",
                "depmap",
                "model_id",
                "cell_line_name",
            )
        ):
            if df[c].nunique(dropna=True) < 100:
                note("MODELS", str(c) + " frequencies", safe_top(df[c], 12))
    columns = [
        c
        for c in df.columns
        if any(
            w in tidy(c)
            for w in (
                "model_id",
                "cosmic",
                "depmap",
                "cell_line_name",
                "sanger",
                "tissue",
                "cancer_type",
                "tcga",
            )
        )
    ]
    if columns:
        save(df[columns].copy(), "model_id_reference.csv")
    result = {}
    for label, aliases in {
        "COSMIC": ["COSMIC_ID", "COSMICID", "COSMIC"],
        "MODEL": ["model_id", "SIDM"],
        "DEPMAP": ["BROAD_ID", "DEPMAP_ID", "ACH_ID", "BROAD"],
        "NAME": ["model_name", "cell_line_name", "cell_line", "model"],
        "TISSUE": ["tissue", "tissue_type", "cancer_type", "tissue_type_detail"],
    }.items():
        col = match_col(df, exact=aliases)
        if col is None and label == "COSMIC":
            col = match_col(df, words=("cosmic",))
        if col is None and label == "DEPMAP":
            col = match_col(df, words=("broad",))
        if col is None and label == "MODEL":
            col = match_col(df, words=("model", "id"))
        if col is None and label == "TISSUE":
            col = match_col(df, words=("tissue",))
        if col is not None:
            result[label] = col
            note("MODELS", label + " field", str(col))
            if label != "TISSUE":
                note("MODELS", label + " unique IDs", clean_id(df[col]).nunique())
    result["_df"] = df
    return result


def gdsc(meta):
    df = read_small("GDSC2_fitted_dose_response_", True)
    if df is None:
        return None
    note("GDSC2", "Rows", f"{len(df):,}")
    note("GDSC2", "Columns", " | ".join(map(str, df.columns)))
    col = {
        "cell": match_col(
            df, exact=["COSMIC_ID", "COSMICID", "CELL_ID", "CELL_LINE_ID"]
        ),
        "drug": match_col(df, exact=["DRUG_ID", "DRUG_NAME", "DRUG"]),
        "name": match_col(df, exact=["DRUG_NAME"]),
        "tissue": match_col(df, exact=["TCGA_DESC", "CANCER_TYPE", "TISSUE"]),
        "ic50": match_col(df, exact=["LN_IC50", "IC50"]),
        "auc": match_col(df, exact=["AUC"]),
        "z": match_col(df, exact=["Z_SCORE", "Z_SCORE_IC50"]),
    }
    for k, v in col.items():
        note("GDSC2", k + " column", str(v))
    if col["cell"]:
        note("GDSC2", "Unique cell lines", f"{clean_id(df[col['cell']]).nunique():,}")
    if col["drug"]:
        note("GDSC2", "Unique drugs", f"{clean_id(df[col['drug']]).nunique():,}")
    if col["tissue"]:
        note("GDSC2", "Cancer types frequency", safe_top(df[col["tissue"]], 18))
    for metric in ("ic50", "auc", "z"):
        if col[metric]:
            note("GDSC2", metric + " distribution", qstats(df[col[metric]]))
            note("GDSC2", metric + " missing", f"{df[col[metric]].isna().sum():,}")
    if col["cell"] and col["drug"]:
        key = (
            clean_id(df[col["cell"]]).fillna("?")
            + "|"
            + clean_id(df[col["drug"]]).fillna("?")
        )
        note("GDSC2", "Repeated cell-drug keys", f"{int(key.duplicated().sum()):,}")
        d = df.copy()
        if col["tissue"] and col["auc"]:
            tmp = (
                d.groupby(col["tissue"], dropna=False)
                .agg(n=(col["auc"], "size"), median_AUC=(col["auc"], "median"))
                .reset_index()
            )
            save(tmp, "gdsc_cancer_response_summary.csv")
        if col["name"] and col["auc"]:
            tmp = (
                d.groupby(col["name"], dropna=False)
                .agg(n=(col["auc"], "size"), median_AUC=(col["auc"], "median"))
                .reset_index()
            )
            tmp = tmp.sort_values("n", ascending=False)
            save(tmp, "gdsc_drug_response_summary.csv")
            note(
                "GDSC2",
                "Most frequently tested drugs",
                "; ".join(str(x) for x in tmp[col["name"]].head(12)),
            )
        if meta and "COSMIC" in meta:
            models = meta["_df"]
            mk = clean_id(models[meta["COSMIC"]])
            dk = clean_id(df[col["cell"]])
            overlap = set(mk.dropna()) & set(dk.dropna())
            note(
                "LINKAGE",
                "GDSC COSMIC IDs matched to model list",
                f"{len(overlap):,} / {dk.nunique():,} ({100 * len(overlap) / max(dk.nunique(), 1):.1f}%)",
            )
            note(
                "LINKAGE",
                "Matched GDSC response rows",
                f"{dk.isin(overlap).sum():,} / {len(df):,}",
            )
    return (df, col)


def prism_one(fragment, section):
    p = getfile(fragment)
    if not p:
        warn(section, "file missing")
        return None
    try:
        head = csv_read(p, nrows=25)
        note(section, "columns", " | ".join(map(str, head.columns)))
        rows = 0
        counts = Counter()
        uniq = defaultdict(set)
        numeric = defaultdict(list)
        id_cols = [
            c
            for c in head.columns
            if any(
                x in tidy(c)
                for x in (
                    "depmap_id",
                    "ccle_name",
                    "cell_line",
                    "drug_id",
                    "compound",
                    "broad_id",
                    "pert_id",
                    "screen_id",
                    "treatment_id",
                )
            )
        ]
        val_cols = [
            c
            for c in head.columns
            if any(
                x in tidy(c)
                for x in (
                    "auc",
                    "ic50",
                    "ec50",
                    "logfold",
                    "log_fold",
                    "viability",
                    "log_ic50",
                    "log_ec50",
                )
            )
        ]
        for chunk in csv_read(p, chunksize=PRISM_CHUNK):
            rows += len(chunk)
            for c in id_cols:
                nonnull = clean_id(chunk[c]).dropna()
                if len(uniq[c]) < 100000:
                    uniq[c].update(nonnull.unique()[:100000])
            for c in val_cols[:8]:
                vals = (
                    pd.to_numeric(chunk[c], errors="coerce")
                    .replace([np.inf, -np.inf], np.nan)
                    .dropna()
                )
                counts[c + "_nonmissing"] += len(vals)
                if len(numeric[c]) < 50000:
                    numeric[c].extend(
                        vals.iloc[: min(len(vals), 50000 - len(numeric[c]))].tolist()
                    )
        note(section, "records", f"{rows:,}")
        for c in id_cols:
            note(
                section,
                "unique " + str(c),
                f"{len(uniq[c]):,}" + (" (capped)" if len(uniq[c]) >= 100000 else ""),
            )
        for c in val_cols[:8]:
            note(
                section,
                str(c) + " [sampled quantiles]",
                qstats(pd.Series(numeric[c]))
                + "; nonmissing="
                + str(counts[c + "_nonmissing"]),
            )
        return dict(header=list(head.columns), ids=uniq, rows=rows)
    except Exception as exc:
        warn(section, exc)
        return None


def prism(gdsc_result, meta):
    prism_one("secondary-screen-dose-response-curve-parameters", "PRISM SECONDARY")
    secondary = read_small("secondary-screen-cell-line-info")
    if secondary is not None:
        note("PRISM META", "Secondary cell info records", len(secondary))
        note(
            "PRISM META",
            "Secondary cell info columns",
            " | ".join(map(str, secondary.columns)),
        )
        for col in secondary.columns:
            if (
                ("depmap" in tidy(col) or "broad" in tidy(col))
                and meta
                and "DEPMAP" in meta
            ):
                mm = clean_id(meta["_df"][meta["DEPMAP"]])
                pp = clean_id(secondary[col])
                both = set(mm.dropna()) & set(pp.dropna())
                note(
                    "LINKAGE",
                    "PRISM cell-line info matched via " + str(col),
                    f"{len(both):,} / {pp.nunique():,} ({100 * len(both) / max(pp.nunique(), 1):.1f}%)",
                )
        if len(secondary) < 100000:
            save(secondary, "prism_secondary_cell_metadata_copy.csv")
    prism_one("primary-screen-replicate-collapsed-logfold-change", "PRISM PRIMARY")
    treatment = read_small("primary-screen-replicate-collapsed-treatment-info")
    if treatment is not None:
        note("PRISM PRIMARY", "Treatment annotation rows", f"{len(treatment):,}")
        note(
            "PRISM PRIMARY",
            "Treatment columns",
            " | ".join(map(str, treatment.columns)),
        )


def omics_summaries(meta):
    specs = [
        ("rnaseq_merged", "RNA EXPRESSION"),
        ("mutations_summary", "DRIVER MUTATIONS"),
        ("cnv_summary", "DRIVER CNV"),
        ("Proteomics_", "PROTEOMICS"),
        ("Project_Score2_fitness_scores", "CRISPR"),
        ("screened_compounds", "DRUG TARGETS"),
        ("gene_identifiers", "GENE ANNOTATION"),
        ("model_dataset_availability", "DATA COVERAGE"),
    ]
    for fragment, section in specs:
        p = getfile(fragment)
        if p is None:
            warn(section, "source not found")
            continue
        try:
            d = peek(p, MAX_INSPECTION_ROWS)
            if d is None:
                note(section, "inspection", "no readable table detected")
                continue
            note(section, "sample rows", len(d))
            note(section, "sample columns", " | ".join(map(str, d.columns[:50])))
            for c in d.columns:
                t = tidy(c)
                if any(
                    w in t
                    for w in (
                        "model_id",
                        "cosmic_id",
                        "depmap_id",
                        "broad_id",
                        "gene_id",
                        "gene_symbol",
                        "tpm",
                        "fpkm",
                        "count",
                        "protein",
                        "alteration",
                        "cancer_type",
                    )
                ):
                    note(
                        section,
                        "sample " + str(c),
                        f"nonnull={d[c].notna().sum()}/{len(d)}; unique={d[c].nunique(dropna=True)}; examples="
                        + ", ".join(map(str, d[c].dropna().astype(str).unique()[:6])),
                    )
        except Exception as exc:
            warn(section, exc)


def main():
    inventory()
    if not DATA.exists():
        return
    meta = model_meta()
    gdsc_result = gdsc(meta)
    prism(gdsc_result, meta)
    omics_summaries(meta)
    note(
        "INTERPRETATION",
        "Experimental response",
        "GDSC LN_IC50 and AUC are continuous drug-response outcomes; they are not classifier sensitivity/specificity.",
    )
    note(
        "INTERPRETATION",
        "Important limitation",
        "Do not compare GDSC2 IC50 with PRISM viability or fitted curve parameters as though they had identical scales.",
    )
    note(
        "INTERPRETATION",
        "Evaluation decision",
        "Define drug/cell/cancer-held-out splits before training, and inspect molecular and drug overlap before claiming external validation.",
    )
    note(
        "INTERPRETATION",
        "Phase 1 scope",
        "Counts/distributions are descriptive; omics ZIP content sampled, NOT comprehensively validated; no predictive model fitted.",
    )


if __name__ == "__main__":
    phase1_failure = None
    try:
        main()
    except Exception as exc:
        phase1_failure = exc
        warn("UNHANDLED FAILURE", traceback.format_exc(limit=6))
    report = pd.DataFrame(REPORT, columns=["section", "metric", "result"])
    save(report, "phase1_report.csv")
    (OUT / "phase1_report.txt").write_text(
        "\n".join(f"[{a}] {b}: {c}" for a, b, c in REPORT), encoding="utf-8"
    )
    print("\n" + "=" * 88)
    print("PROJECT 5 | SOURCE DATA AUDIT")
    print("=== BEGIN PHASE 1 RESULTS ===")
    for a, b, c in REPORT:
        print(f"[{a}] {b}: {c}")
    print(f"\n[COMPLETION] Issues/warnings: {len(ISSUES)}")
    print("[COMPLETION] Output folder:", OUT)
    print("=== END PHASE 1 RESULTS ===")
    print("=" * 88)
    if phase1_failure is not None:
        raise phase1_failure

# ----------------------------------------------------------------------------
# Multi-omics availability and identifier harmonisation
# ----------------------------------------------------------------------------

# ============================================================
# PROJECT 5 - PHASE 1B
# DATASET INTEGRATION AND MODEL READINESS
#
# Python 3.13+
# Required: pandas, numpy, openpyxl
#
# Input:
# Project5/Data
#
# Output:
# Project5/Results/Phase1
#
# All console results are printed at the END.
# ============================================================

from pathlib import Path
from collections import Counter
from contextlib import contextmanager
import zipfile
import re
import traceback

import numpy as np
import pandas as pd


# ============================================================
# 1. CONFIGURATION
# ============================================================

ROOT = PROJECT_ROOT
DATA = ROOT / "Data"
OUT = ROOT / "Results" / "Phase1"

OUT.mkdir(parents=True, exist_ok=True)

REPORT = []
ERRORS = []


def log(section, message):
    REPORT.append(f"[{section}] {message}")


def report_error(section, error):
    ERRORS.append(f"{section}: {error}")
    log("WARNING", f"{section}: {error}")


def find_file(prefix):
    matches = sorted(
        p
        for p in DATA.iterdir()
        if p.is_file() and p.name.lower().startswith(prefix.lower())
    )

    if not matches:
        raise FileNotFoundError(f"Missing file beginning with: {prefix}")

    return matches[0]


def ids(series):
    return (
        series.astype("string")
        .str.strip()
        .str.upper()
        .str.replace(r"\.0$", "", regex=True)
        .replace("", pd.NA)
        .dropna()
    )


def unique_ids(series):
    return set(ids(series))


def percentage(n, total):
    if total == 0:
        return "N/A"
    return f"{100 * n / total:.2f}%"


def normalize_name(series):
    return (
        series.astype("string")
        .str.upper()
        .str.strip()
        .str.replace(r"[^A-Z0-9]", "", regex=True)
        .replace("", pd.NA)
    )


def require_columns(df, columns, dataset):
    missing = [c for c in columns if c not in df.columns]

    if missing:
        raise ValueError(f"{dataset}: missing columns {missing}")


def save(df, filename):
    df.to_csv(OUT / filename, index=False)


# ============================================================
# 2. READ COMPRESSED DATA WITHOUT EXTRACTING ZIP FILES
# ============================================================


@contextmanager
def zip_csv_stream(zip_path, preferred=None):
    with zipfile.ZipFile(zip_path) as archive:
        members = [
            name
            for name in archive.namelist()
            if name.lower().endswith((".csv", ".tsv"))
        ]

        if preferred:
            selected = [name for name in members if preferred.lower() in name.lower()]

            if not selected:
                raise ValueError(
                    f"Member containing '{preferred}' not found in {zip_path.name}"
                )

            members = selected

        if not members:
            raise ValueError(f"No CSV/TSV found inside {zip_path.name}")

        # Prefer smaller files unless a specific member was requested.
        member = min(members, key=lambda x: archive.getinfo(x).file_size)

        log("ZIP", f"{zip_path.name}: selected {member}")

        with archive.open(member) as stream:
            yield stream, member


def read_zip_csv(zip_path, preferred=None, **kwargs):
    with zip_csv_stream(zip_path, preferred) as (stream, member):
        sep = "\t" if member.lower().endswith(".tsv") else ","

        return pd.read_csv(stream, sep=sep, low_memory=False, **kwargs)


def zip_header(zip_path, preferred):
    return read_zip_csv(zip_path, preferred=preferred, nrows=0).columns.tolist()


# ============================================================
# 3. LOAD MODEL METADATA
# ============================================================

log("START", "Loading Cell Model Passports metadata")

models = pd.read_csv(find_file("model_list_"), low_memory=False)

require_columns(
    models, ["model_id", "model_name", "tissue", "cancer_type"], "Model metadata"
)

models["model_id"] = ids(models["model_id"]).reindex(models.index)

models = models.dropna(subset=["model_id"]).copy()

if models["model_id"].duplicated().any():
    raise ValueError("Duplicate model IDs detected in model list")

model_ids = set(models["model_id"])

log("MODELS", f"Total models: {len(models):,}")
log("MODELS", f"Unique model IDs: {len(model_ids):,}")

if "BROAD_ID" in models:
    broad = ids(models["BROAD_ID"])

    log("MODELS", f"Models with Broad identifiers: {len(broad):,}")

log(
    "MODELS",
    "Top cancer types: " + str(models["cancer_type"].value_counts().head(12).to_dict()),
)


# ============================================================
# 4. GDSC2 IDENTIFIER CORRECTION
# ============================================================

log("GDSC2", "Loading fitted drug-response data")

gdsc = pd.read_excel(find_file("GDSC2_fitted_dose_response_"))

require_columns(
    gdsc, ["SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "LN_IC50", "AUC"], "GDSC2"
)

gdsc["model_id"] = ids(gdsc["SANGER_MODEL_ID"]).reindex(gdsc.index)

gdsc["DRUG_ID"] = ids(gdsc["DRUG_ID"]).reindex(gdsc.index)

gdsc_models = unique_ids(gdsc["model_id"])

matched_gdsc = gdsc_models & model_ids
unmatched_gdsc = gdsc_models - model_ids

log("GDSC2", f"Total response records: {len(gdsc):,}")

log("GDSC2", f"Unique cell lines: {len(gdsc_models):,}")

log("GDSC2", f"Unique drugs: {gdsc['DRUG_ID'].nunique():,}")

log("GDSC2", f"Matched Sanger models: {len(matched_gdsc):,}")

log("GDSC2", f"Matching percentage: {percentage(len(matched_gdsc), len(gdsc_models))}")

log("GDSC2", f"Unmatched models: {len(unmatched_gdsc):,}")

gdsc_matched = gdsc[gdsc["model_id"].isin(matched_gdsc)].copy()

gdsc_matched = gdsc_matched.merge(
    models[["model_id", "model_name", "tissue", "cancer_type"]],
    on="model_id",
    how="left",
    validate="many_to_one",
    suffixes=("", "_MODEL"),
)

log("GDSC2", f"Matched response records: {len(gdsc_matched):,}")

# Inspect technical replication: do not silently treat repeated
# model-drug measurements as independent samples.

pair_counts = gdsc_matched.groupby(["model_id", "DRUG_ID"]).size()

log("GDSC2", f"Unique model-drug pairs: {len(pair_counts):,}")

log("GDSC2", f"Pairs with multiple records: {int((pair_counts > 1).sum()):,}")

log("GDSC2", f"Maximum records per pair: {pair_counts.max()}")

save(
    pd.DataFrame(
        {
            "model_id": sorted(gdsc_models),
            "matched": [x in matched_gdsc for x in sorted(gdsc_models)],
        }
    ),
    "gdsc_model_matching.csv",
)

del gdsc_matched


# ============================================================
# 5. RNA-SEQ EXPRESSION COVERAGE
# ============================================================

log("RNA", "Inspecting actual expression matrix")

rna_zip = find_file("rnaseq_merged_")

# Use the wide RSEM TPM matrix, not the 5 GB long-format file.
rna_file = "rnaseq_merged_rsem_tpm_"

rna_sample = read_zip_csv(rna_zip, preferred=rna_file, nrows=8)

rna_columns = rna_sample.columns.tolist()

rna_model_columns = {
    str(c).strip().upper()
    for c in rna_columns
    if re.fullmatch(r"SIDM\d+", str(c).strip(), re.I)
}

log("RNA", f"Expression matrix columns: {len(rna_columns):,}")

log("RNA", f"Cell-line columns: {len(rna_model_columns):,}")

log(
    "RNA",
    "First-column preview: " + str(rna_sample.iloc[:8, :4].to_dict("records"))[:1200],
)

rna_gdsc = rna_model_columns & gdsc_models

log("RNA", f"RNA-seq + GDSC2 shared models: {len(rna_gdsc):,}")

log("RNA", f"GDSC2 models with RNA-seq: {percentage(len(rna_gdsc), len(gdsc_models))}")

log(
    "RNA",
    "Note: Expression matrix orientation and metadata rows "
    "must be handled before feature extraction.",
)


# ============================================================
# 6. MUTATION COVERAGE
# ============================================================

mutations = read_zip_csv(
    find_file("mutations_summary_"), preferred="mutations_summary_"
)

require_columns(mutations, ["model_id", "gene_symbol"], "Mutation summary")

mutation_models = unique_ids(mutations["model_id"])

log("MUTATIONS", f"Mutation records: {len(mutations):,}")

log("MUTATIONS", f"Models with mutation records: {len(mutation_models):,}")

log("MUTATIONS", f"Shared with GDSC2: {len(mutation_models & gdsc_models):,}")

log(
    "MUTATIONS",
    "Important: absence from this mutation table "
    "does not automatically establish wild-type status.",
)


# ============================================================
# 7. CNV COVERAGE
# ============================================================

cnv_zip = find_file("cnv_summary_")

cnv_models = set()
cnv_records = 0

with zip_csv_stream(cnv_zip, "cnv_summary_") as (stream, member):
    for chunk in pd.read_csv(stream, usecols=["model_id"], chunksize=250_000):
        cnv_records += len(chunk)
        cnv_models.update(unique_ids(chunk["model_id"]))

log("CNV", f"CNV records: {cnv_records:,}")

log("CNV", f"Models represented: {len(cnv_models):,}")

log("CNV", f"Shared with GDSC2: {len(cnv_models & gdsc_models):,}")


# ============================================================
# 8. PROTEOMICS COVERAGE
# ============================================================

protein_zip = find_file("Proteomics_")

protein_models = set()
protein_records = 0

with zip_csv_stream(protein_zip, "proteomics_all_") as (stream, member):
    for chunk in pd.read_csv(stream, usecols=["model_id"], chunksize=250_000):
        protein_records += len(chunk)
        protein_models.update(unique_ids(chunk["model_id"]))

log("PROTEOMICS", f"Protein records: {protein_records:,}")

log("PROTEOMICS", f"Models represented: {len(protein_models):,}")

log("PROTEOMICS", f"Shared with GDSC2: {len(protein_models & gdsc_models):,}")


# ============================================================
# 9. CRISPR COVERAGE
# ============================================================

crispr_zip = find_file("Project_Score2_")

crispr_members = []

with zipfile.ZipFile(crispr_zip) as z:
    crispr_members = [name for name in z.namelist() if name.endswith(".tsv")]

log("CRISPR", f"Available matrix files: {len(crispr_members)}")

for member in crispr_members:
    header = read_zip_csv(crispr_zip, preferred=member.split("/")[-1], nrows=0)

    log("CRISPR", f"{Path(member).name}: {len(header.columns):,} columns")

log(
    "CRISPR",
    "Model-name mapping is required before CRISPR "
    "coverage can be considered confirmed.",
)


# ============================================================
# 10. MULTI-OMICS OVERLAP
# ============================================================

coverage = pd.DataFrame({"model_id": sorted(gdsc_models)})

coverage["RNA"] = coverage["model_id"].isin(rna_model_columns)
coverage["MUTATION_RECORD"] = coverage["model_id"].isin(mutation_models)
coverage["CNV"] = coverage["model_id"].isin(cnv_models)
coverage["PROTEOMICS"] = coverage["model_id"].isin(protein_models)

log("COVERAGE", f"GDSC2 models: {len(coverage):,}")

for col in ["RNA", "MUTATION_RECORD", "CNV", "PROTEOMICS"]:
    n = int(coverage[col].sum())

    log("COVERAGE", f"{col}: {n:,} ({percentage(n, len(coverage))})")

for combo in [["RNA", "CNV"], ["RNA", "PROTEOMICS"], ["RNA", "CNV", "PROTEOMICS"]]:
    n = int(coverage[combo].all(axis=1).sum())

    log("COVERAGE", f"{' + '.join(combo)}: {n:,}")

save(coverage, "gdsc_multiomics_coverage.csv")


# ============================================================
# 11. MODEL AVAILABILITY METADATA
# ============================================================

availability = pd.read_csv(find_file("model_dataset_availability_"), low_memory=False)

require_columns(availability, ["model_id"], "Dataset availability")

availability["model_id"] = ids(availability["model_id"]).reindex(availability.index)

availability = availability[availability["model_id"].isin(gdsc_models)]

log(
    "AVAILABILITY",
    f"Models represented in availability metadata: "
    f"{availability['model_id'].nunique():,}",
)

for column in availability.columns:
    if any(
        term in column.lower()
        for term in ["rnaseq", "mutation", "proteomics", "crispr"]
    ):
        counts = availability[column].value_counts(dropna=False).head(5).to_dict()

        log("AVAILABILITY", f"{column}: {counts}")


# ============================================================
# 12. PRISM CELL-LINE MATCHING
# ============================================================

prism = pd.read_csv(
    find_file("secondary-screen-dose-response-curve-parameters"),
    usecols=["broad_id", "depmap_id", "name", "auc", "ic50", "screen_id"],
    low_memory=False,
)

prism["depmap_id"] = ids(prism["depmap_id"]).reindex(prism.index)

prism_models = unique_ids(prism["depmap_id"])

if "BROAD_ID" in models:
    model_bridge = models[["model_id", "BROAD_ID"]].dropna().copy()

    model_bridge["BROAD_ID"] = ids(model_bridge["BROAD_ID"]).reindex(model_bridge.index)

    model_bridge = model_bridge.dropna().drop_duplicates()

    ambiguous = model_bridge.groupby("BROAD_ID")["model_id"].nunique()

    ambiguous_ids = set(ambiguous[ambiguous > 1].index)

    if ambiguous_ids:
        log("PRISM", f"Ambiguous Broad IDs excluded: {len(ambiguous_ids):,}")
        model_bridge = model_bridge[~model_bridge["BROAD_ID"].isin(ambiguous_ids)]

    bridge_lookup = model_bridge.drop_duplicates("BROAD_ID").set_index("BROAD_ID")[
        "model_id"
    ]

    shared_prism = prism_models & set(bridge_lookup.index)

    matched_prism_sanger = {bridge_lookup[x] for x in shared_prism}

    cross_platform_cells = matched_prism_sanger & gdsc_models

    log("PRISM", f"Unique PRISM models: {len(prism_models):,}")

    log("PRISM", f"PRISM models matched to Sanger: {len(shared_prism):,}")

    log("PRISM", f"Cell lines shared with GDSC2: {len(cross_platform_cells):,}")

    save(model_bridge, "sanger_prism_model_bridge.csv")

else:
    cross_platform_cells = set()
    log("PRISM", "BROAD_ID mapping unavailable")


# ============================================================
# 13. DRUG MATCHING - GDSC2 VERSUS PRISM
# ============================================================

gdsc_drugs = gdsc[["DRUG_ID", "DRUG_NAME"]].drop_duplicates().copy()

prism_drugs = prism[["broad_id", "name"]].drop_duplicates().copy()

gdsc_drugs["match_name"] = normalize_name(gdsc_drugs["DRUG_NAME"])

prism_drugs["match_name"] = normalize_name(prism_drugs["name"])

gdsc_named = gdsc_drugs.dropna(subset=["match_name"])
prism_named = prism_drugs.dropna(subset=["match_name"])

# Conservative first-pass matching by normalised names.
# These matches are candidates, not yet validated chemical identities.

drug_matches = gdsc_named.merge(prism_named, on="match_name", how="inner")

log("DRUG MATCHING", f"GDSC2 unique drug IDs: {gdsc['DRUG_ID'].nunique():,}")

log("DRUG MATCHING", f"PRISM unique compound IDs: {prism['broad_id'].nunique():,}")

log("DRUG MATCHING", f"Candidate matched compound pairs: {len(drug_matches):,}")

log(
    "DRUG MATCHING",
    f"Unique GDSC2 drugs with name matches: {drug_matches['DRUG_ID'].nunique():,}",
)

log(
    "DRUG MATCHING",
    "Matches are provisional: chemical identity and "
    "formulation must be verified before external validation.",
)

save(drug_matches, "gdsc_prism_candidate_drug_matches.csv")


# ============================================================
# 14. JOINT MATCHING COVERAGE
# ============================================================

if cross_platform_cells:
    joint = cross_platform_cells & rna_model_columns

    log("JOINT COVERAGE", f"GDSC2 + PRISM + RNA models: {len(joint):,}")

    joint_protein = joint & protein_models

    log("JOINT COVERAGE", f"GDSC2 + PRISM + RNA + proteomics: {len(joint_protein):,}")


# ============================================================
# 15. MODELLING READINESS
# ============================================================

log("DECISION", "Primary outcome: continuous drug response.")

log("DECISION", "Primary input: gene expression.")

log(
    "DECISION",
    "Secondary inputs: mutations, CNV and proteomics "
    "where sufficient matched coverage exists.",
)

log(
    "DECISION",
    "Train/validation/test groups must be assigned "
    "before model fitting and feature selection.",
)

log(
    "DECISION",
    "Evaluate held-out cell lines and cancer types; "
    "consider held-out drugs only when suitable drug "
    "descriptors are available.",
)

log(
    "DECISION",
    "For drug-response classification, thresholds must "
    "be scientifically justified and defined using "
    "training data only.",
)

log(
    "DECISION",
    "PRISM is a separate experimental platform; "
    "its response values must not be treated as "
    "equivalent to GDSC2 outcomes.",
)

log(
    "DECISION",
    "Graph models will be considered after "
    "establishing conventional prediction baselines.",
)


# ============================================================
# 16. SAVE RESULTS
# ============================================================

report_path = OUT / "phase1_integration_report.txt"

with open(report_path, "w", encoding="utf-8") as f:
    f.write("\n".join(REPORT))

log("COMPLETION", f"Results folder: {OUT}")


# ============================================================
# 17. PRINT ALL OUTPUTS AT THE END
# ============================================================

print("\n" + "=" * 90)

print("PROJECT 5 - PHASE 1B | DATA INTEGRATION RESULTS")

print("=" * 90)

print("\n=== BEGIN PHASE 1B RESULTS ===\n")

for line in REPORT:
    print(line)

print("\n" + "-" * 90)

print(f"TOTAL WARNINGS: {len(ERRORS)}")

if ERRORS:
    for error in ERRORS:
        print(error)

print("\n=== END PHASE 1B RESULTS ===")

print("=" * 90)

if ERRORS:
    raise RuntimeError(f"Phase 1B completed with {len(ERRORS)} error(s); review the printed report")
