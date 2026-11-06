"""Validate the six source specimens and prepare harmonised count matrices.

The script audits the supplied archives, checks barcode and metadata alignment,
applies the established exclusions and writes analysis-ready sparse matrices.
Place the three source archives in a ``Data`` directory beside this file.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ============================================================================
# Inspect archive contents and source-file structure
# ============================================================================

from pathlib import Path
import tarfile
import csv
import io
import re
from collections import Counter, defaultdict

# ============================================================
# PROJECT 2 | STAGE 1A
# Wu et al. (2021) — Breast Cancer Spatial Transcriptomics
# Dataset inventory and initial inspection
# ============================================================


ROOT = PROJECT_ROOT / "Data"

OUTPUT = ROOT.parent / "Results" / "Stage1A"
OUTPUT.mkdir(parents=True, exist_ok=True)

ARCHIVES = ["filtered_count_matrices.tar.gz", "spatial.tar.gz", "metadata.tar.gz"]

MAX_PREVIEW_ROWS = 5
MAX_DISPLAY_MEMBERS = 30

report = []


def log(message=""):
    report.append(str(message))


def format_size(size):
    return f"{size / 1024**2:,.2f} MB"


def inspect_table(archive, member):
    """
    Preview small tabular files without extracting archives.
    Large expression files are not loaded into memory.
    """

    try:
        stream = archive.extractfile(member)

        if stream is None:
            return

        # Only preview the beginning of each file.
        content = stream.read(65536).decode("utf-8-sig", errors="replace")

        if member.name.lower().endswith(".csv"):
            delimiter = ","
        else:
            delimiter = "\t"

        reader = csv.reader(io.StringIO(content), delimiter=delimiter)

        rows = []

        for i, row in enumerate(reader):
            rows.append(row)

            if i >= MAX_PREVIEW_ROWS:
                break

        log(f"\nTABLE PREVIEW: {member.name}")

        if rows:
            log(f"Columns detected: {len(rows[0])}")
            log(f"First row: {rows[0][:15]}")

            for row in rows[1:]:
                log(f"  {row[:15]}")

    except Exception as error:
        log(f"Preview error: {error}")


def inspect_matrix_header(archive, member):
    """
    Read Matrix Market dimensions without loading counts.
    """

    try:
        stream = archive.extractfile(member)

        if stream is None:
            return

        for line in stream:
            decoded = line.decode("utf-8", errors="replace").strip()

            if not decoded or decoded.startswith("%"):
                continue

            values = decoded.split()

            if len(values) >= 3:
                log(
                    f"Matrix dimensions: "
                    f"{values[0]} rows × "
                    f"{values[1]} columns; "
                    f"{values[2]} stored entries"
                )

            break

    except Exception as error:
        log(f"Matrix header error: {error}")


# ============================================================
# 1. CHECK INPUT DIRECTORY
# ============================================================

log("=" * 75)
log("PROJECT 2 | STAGE 1A — DATASET INSPECTION")
log("=" * 75)

log(f"\nInput directory: {ROOT}")
log(f"Directory exists: {ROOT.exists()}")

if not ROOT.exists():
    raise FileNotFoundError(ROOT)

log("\nFILES IN DATA DIRECTORY")

for path in sorted(ROOT.iterdir()):
    if path.is_file():
        log(f"{path.name} | {format_size(path.stat().st_size)}")

    elif path.is_dir():
        log(f"{path.name} | DIRECTORY")


# ============================================================
# 2. INSPECT ARCHIVES
# ============================================================

archive_summaries = []

for archive_name in ARCHIVES:
    path = ROOT / archive_name

    log("\n" + "=" * 75)
    log(f"ARCHIVE: {archive_name}")
    log("=" * 75)

    if not path.exists():
        log("STATUS: NOT FOUND")

        archive_summaries.append(
            {"archive": archive_name, "status": "MISSING", "files": 0, "directories": 0}
        )

        continue

    try:
        with tarfile.open(path, "r:gz") as tar:
            members = tar.getmembers()

            files = [m for m in members if m.isfile()]
            directories = [m for m in members if m.isdir()]

            log("STATUS: READABLE")
            log(f"Files: {len(files):,}")
            log(f"Directories: {len(directories):,}")

            total_size = sum(m.size for m in files)

            log(f"Uncompressed file size: {format_size(total_size)}")

            archive_summaries.append(
                {
                    "archive": archive_name,
                    "status": "READABLE",
                    "files": len(files),
                    "directories": len(directories),
                }
            )

            # ----------------------------------------
            # File extensions
            # ----------------------------------------

            extensions = Counter()

            for member in files:
                name = member.name.lower()

                if name.endswith(".mtx.gz"):
                    ext = ".mtx.gz"
                elif name.endswith(".tsv.gz"):
                    ext = ".tsv.gz"
                elif name.endswith(".csv.gz"):
                    ext = ".csv.gz"
                else:
                    ext = Path(name).suffix or "[none]"

                extensions[ext] += 1

            log("\nFILE TYPES")

            for ext, count in extensions.most_common():
                log(f"{ext}: {count}")

            # ----------------------------------------
            # Directory structure
            # ----------------------------------------

            log("\nFIRST ARCHIVE MEMBERS")

            for member in members[:MAX_DISPLAY_MEMBERS]:
                kind = "DIR " if member.isdir() else "FILE"

                log(f"[{kind}] {member.name}")

            if len(members) > MAX_DISPLAY_MEMBERS:
                log(f"... {len(members) - MAX_DISPLAY_MEMBERS} additional members")

            # ----------------------------------------
            # Potential sample directories
            # ----------------------------------------

            top_folders = Counter()

            for member in files:
                parts = Path(member.name).parts

                if len(parts) >= 2:
                    top_folders["/".join(parts[:2])] += 1

            log("\nTOP DIRECTORY GROUPS")

            for folder, count in top_folders.most_common(20):
                log(f"{folder}: {count} files")

            # ----------------------------------------
            # Identify important files
            # ----------------------------------------

            keywords = [
                "matrix",
                "barcode",
                "feature",
                "gene",
                "tissue_positions",
                "scalefactor",
                "metadata",
                "annotation",
                "pathology",
            ]

            important = [
                m for m in files if any(key in m.name.lower() for key in keywords)
            ]

            log("\nIMPORTANT FILES")

            for member in important[:40]:
                log(member.name)

            if len(important) > 40:
                log(f"... {len(important) - 40} additional matches")

            # ----------------------------------------
            # Matrix Market dimensions
            # ----------------------------------------

            matrix_files = [m for m in files if m.name.lower().endswith(".mtx")]

            log("\nMATRIX MARKET FILES")

            for member in matrix_files[:20]:
                log(f"\n{member.name}")

                inspect_matrix_header(tar, member)

            # ----------------------------------------
            # Small metadata tables
            # ----------------------------------------

            table_files = [
                m
                for m in files
                if m.name.lower().endswith((".csv", ".tsv")) and m.size <= 5 * 1024**2
            ]

            log("\nSMALL TABULAR FILES")

            for member in table_files[:15]:
                log(f"{member.name} | {format_size(member.size)}")

                inspect_table(tar, member)

    except Exception as error:
        log(f"STATUS: ERROR — {error}")

        archive_summaries.append(
            {
                "archive": archive_name,
                "status": f"ERROR: {error}",
                "files": 0,
                "directories": 0,
            }
        )


# ============================================================
# 3. FINAL CONSOLIDATED SUMMARY
# ============================================================

log("\n" + "=" * 75)
log("FINAL SUMMARY")
log("=" * 75)

for item in archive_summaries:
    log(
        f"{item['archive']} | "
        f"{item['status']} | "
        f"{item['files']:,} files | "
        f"{item['directories']:,} directories"
    )

readable = sum(item["status"] == "READABLE" for item in archive_summaries)

log(f"\nReadable archives: {readable}/{len(ARCHIVES)}")

log("\nSTAGE 1A CONCLUSION")

if readable == len(ARCHIVES):
    log("All required archives are accessible.")
    log("Next: inspect actual expression and spatial schemas.")
else:
    log("Some archives are missing or unreadable.")
    log("Resolve these issues before advancing.")


# ============================================================
# 4. SAVE AND PRINT COMPLETE REPORT
# ============================================================

report_path = OUTPUT / "Stage1A_Inspection_Report.txt"

report_path.write_text("\n".join(report), encoding="utf-8")

print("\n".join(report))

print("\n" + "=" * 75)
print("REPORT SAVED")
print("=" * 75)

print(report_path)

print("\nSTAGE 1A COMPLETE")

# ============================================================================
# Validate expression matrices, coordinates and pathology metadata
# ============================================================================

from pathlib import Path
from collections import Counter
import gzip
import io
import json
import tarfile
import traceback

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse


ROOT = PROJECT_ROOT / "Data"
OUT = ROOT.parent / "Results" / "Stage1B"
OUT.mkdir(parents=True, exist_ok=True)

COUNT_ARCHIVE = ROOT / "filtered_count_matrices.tar.gz"
SPACE_ARCHIVE = ROOT / "spatial.tar.gz"
META_ARCHIVE = ROOT / "metadata.tar.gz"
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

summary = []
checks = []
spot_qc = []
annotation_counts = []
errors = []


def say(s=""):
    summary.append(str(s))


def get_member(tar, name):
    try:
        member = tar.getmember(name)
    except KeyError as exc:
        raise FileNotFoundError(f"Missing archive member: {name}") from exc
    handle = tar.extractfile(member)
    if handle is None:
        raise ValueError(f"Cannot read archive member: {name}")
    return handle


def read_maybe_gzip(tar, name):
    """Open archive member; detect gzip by magic bytes, not .gz suffix.

    Wu archive contains some files labelled .gz whose bytes are plain text.
    """
    source = get_member(tar, name)
    wrapped = io.BufferedReader(source)
    if wrapped.peek(2)[:2] == b"\x1f\x8b":
        return gzip.GzipFile(fileobj=wrapped)
    return wrapped


def tsv_gz(tar, name, columns):
    with read_maybe_gzip(tar, name) as stream:
        data = pd.read_csv(stream, sep="\t", header=None, dtype=str)
    if data.shape[1] > len(columns):
        data = data.iloc[:, : len(columns)]
    if data.shape[1] == 0:
        raise ValueError(f"{name}: empty feature/barcode table")
    data.columns = columns[: data.shape[1]]
    for col in columns[data.shape[1] :]:
        data[col] = pd.NA  # e.g. older 10x features.tsv has no feature_type
    return data[columns]


def percentile(values, q):
    return float(np.percentile(values, q)) if len(values) else float("nan")


say("PROJECT 2 | STAGE 1B — EXPRESSION, SPATIAL & PATHOLOGY VALIDATION")
say("Wu et al. (2021) | Six breast cancer Visium specimens")
say("No QC thresholds applied; no spots removed; raw counts preserved.")
say("")

for item in (COUNT_ARCHIVE, SPACE_ARCHIVE, META_ARCHIVE):
    if not item.is_file():
        raise FileNotFoundError(f"Required archive not found: {item}")

with (
    tarfile.open(COUNT_ARCHIVE, "r:gz") as counts_tar,
    tarfile.open(SPACE_ARCHIVE, "r:gz") as spatial_tar,
    tarfile.open(META_ARCHIVE, "r:gz") as meta_tar,
):
    for sample in SAMPLES:
        say("=" * 70)
        say(f"SAMPLE: {sample}")
        try:
            prefix = f"filtered_count_matrices/{sample}_filtered_count_matrix"
            feat = tsv_gz(
                counts_tar,
                f"{prefix}/features.tsv.gz",
                ["gene_id", "gene_symbol", "feature_type"],
            )
            bar = tsv_gz(counts_tar, f"{prefix}/barcodes.tsv.gz", ["barcode"])
            genes = feat["gene_symbol"].fillna("").astype(str)
            barcode = bar["barcode"].astype(str)

            with read_maybe_gzip(counts_tar, f"{prefix}/matrix.mtx.gz") as stream:
                matrix = mmread(stream)
            matrix = (
                matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
            )
            matrix.sum_duplicates()

            if matrix.shape != (len(feat), len(bar)):
                raise ValueError(
                    f"Matrix shape {matrix.shape} incompatible with features {len(feat)} and barcodes {len(bar)}"
                )
            if barcode.duplicated().any():
                raise ValueError("Duplicate expression barcodes")
            if (matrix.data < 0).any():
                raise ValueError("Negative counts detected")

            spatial_name = f"spatial/{sample}_spatial/tissue_positions_list.csv"
            with get_member(spatial_tar, spatial_name) as src:
                positions = pd.read_csv(
                    src,
                    header=None,
                    names=[
                        "barcode",
                        "in_tissue",
                        "array_row",
                        "array_col",
                        "pixel_row",
                        "pixel_col",
                    ],
                    dtype={"barcode": str},
                )
            if positions["barcode"].duplicated().any():
                raise ValueError("Duplicate spatial barcodes")
            positions["in_tissue"] = pd.to_numeric(
                positions["in_tissue"], errors="coerce"
            )

            with get_member(
                spatial_tar, f"spatial/{sample}_spatial/scalefactors_json.json"
            ) as src:
                scalefactors = json.load(src)

            with get_member(meta_tar, f"metadata/{sample}_metadata.csv") as src:
                meta = pd.read_csv(src, index_col=0)
            meta.index = meta.index.astype(str)
            if meta.index.duplicated().any():
                raise ValueError("Duplicate metadata barcodes")
            mandatory = [
                "nCount_RNA",
                "nFeature_RNA",
                "patientid",
                "subtype",
                "Classification",
            ]
            missing_cols = set(mandatory).difference(meta.columns)
            if missing_cols:
                raise ValueError(f"Missing metadata columns: {sorted(missing_cols)}")

            # Matrix is genes x spots; matrix operations return per-spot vectors.
            total_counts = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
            genes_detected = np.asarray(matrix.getnnz(axis=0)).ravel().astype(np.int64)
            qc = pd.DataFrame(
                {
                    "sample": sample,
                    "barcode": barcode.to_numpy(),
                    "total_counts": total_counts,
                    "genes_detected": genes_detected,
                }
            )
            coords = positions.set_index("barcode")
            aligned = qc.join(coords, on="barcode")
            aligned = aligned.join(meta, on="barcode")
            expr_barcodes = set(barcode)
            pos_barcodes = set(positions["barcode"])
            meta_barcodes = set(meta.index)
            common = expr_barcodes & pos_barcodes & meta_barcodes
            pos_missing = int((~qc["barcode"].isin(pos_barcodes)).sum())
            meta_missing = int((~qc["barcode"].isin(meta_barcodes)).sum())
            metadata_total = pd.to_numeric(aligned["nCount_RNA"], errors="coerce")
            metadata_genes = pd.to_numeric(aligned["nFeature_RNA"], errors="coerce")
            count_mismatch = int(
                ((metadata_total - aligned["total_counts"]).abs() > 0.5).sum()
            )
            gene_mismatch = int(
                ((metadata_genes - aligned["genes_detected"]).abs() > 0.5).sum()
            )
            subtype_values = sorted(
                meta["subtype"].dropna().astype(str).unique().tolist()
            )
            patient_values = sorted(
                meta["patientid"].dropna().astype(str).unique().tolist()
            )
            class_counts = (
                meta["Classification"]
                .fillna("[MISSING]")
                .astype(str)
                .value_counts(dropna=False)
            )

            annotation_counts.extend(
                [
                    {"sample": sample, "classification": key, "n_spots": int(value)}
                    for key, value in class_counts.items()
                ]
            )
            spot_qc.append(aligned)
            checks.append(
                {
                    "sample": sample,
                    "genes": matrix.shape[0],
                    "expression_spots": matrix.shape[1],
                    "spatial_spots": len(positions),
                    "metadata_spots": len(meta),
                    "all_three_matched": len(common),
                    "missing_spatial_from_expression": pos_missing,
                    "missing_metadata_from_expression": meta_missing,
                    "nCount_mismatches": count_mismatch,
                    "nFeature_mismatches": gene_mismatch,
                    "gene_id_duplicates": int(feat["gene_id"].duplicated().sum()),
                    "symbol_duplicates": int(genes[genes.ne("")].duplicated().sum()),
                    "raw_nonzero_entries": int(matrix.nnz),
                    "tissue_spots_in_matrix": int((aligned["in_tissue"] == 1).sum()),
                    "subtypes": ";".join(subtype_values),
                    "patients": ";".join(patient_values),
                }
            )

            say(
                f"Expression: {matrix.shape[0]:,} genes × {matrix.shape[1]:,} spots; nonzero entries {matrix.nnz:,}"
            )
            say(
                f"Coordinates: {len(positions):,} spots; metadata: {len(meta):,}; all-three matched: {len(common):,}"
            )
            say(
                f"Tissue spots in expression matrix: {(aligned['in_tissue'] == 1).sum():,}"
            )
            say(
                f"Median counts/spot: {percentile(total_counts, 50):,.0f} [Q1 {percentile(total_counts, 25):,.0f}, Q3 {percentile(total_counts, 75):,.0f}]"
            )
            say(
                f"Median detected genes/spot: {percentile(genes_detected, 50):,.0f} [Q1 {percentile(genes_detected, 25):,.0f}, Q3 {percentile(genes_detected, 75):,.0f}]"
            )
            say(f"Patient IDs: {patient_values}; subtypes: {subtype_values}")
            say(
                f"Metadata count/gene discrepancies: {count_mismatch}/{gene_mismatch} (counts compared against full raw matrix)"
            )
            say(
                f"Expression spots missing coordinates/metadata: {pos_missing}/{meta_missing}"
            )
            say(f"Spatial scale factors: {scalefactors}")
            say("Pathology annotations:")
            for key, value in class_counts.items():
                say(f"  {key}: {value:,}")
        except Exception as exc:
            errors.append({"sample": sample, "error": f"{type(exc).__name__}: {exc}"})
            say(f"ERROR: {type(exc).__name__}: {exc}")
            say(traceback.format_exc())

if checks:
    pd.DataFrame(checks).to_csv(OUT / "sample_integrity_summary.csv", index=False)
if annotation_counts:
    pd.DataFrame(annotation_counts).to_csv(
        OUT / "pathology_annotation_counts.csv", index=False
    )
if spot_qc:
    pd.concat(spot_qc, ignore_index=True).to_csv(
        OUT / "spot_level_qc_unfiltered.csv", index=False
    )
if errors:
    pd.DataFrame(errors).to_csv(OUT / "errors.csv", index=False)

say("=" * 70)
say("FINAL CONSOLIDATED RESULTS")
say(f"Successfully inspected: {len(checks)}/{len(SAMPLES)} samples")
if checks:
    check_df = pd.DataFrame(checks)
    say(f"Total expression spots: {check_df['expression_spots'].sum():,}")
    say(
        f"Total expression spots with both spatial and metadata records: {check_df['all_three_matched'].sum():,}"
    )
    say(
        f"Distinct gene counts across samples: {sorted(check_df['genes'].unique().tolist())}"
    )
    say(f"Total metadata count discrepancies: {check_df['nCount_mismatches'].sum():,}")
    say(
        f"Total metadata feature discrepancies: {check_df['nFeature_mismatches'].sum():,}"
    )
    say("Sample integrity table:")
    say(check_df.to_string(index=False))
if errors:
    say("ERRORS:")
    for err in errors:
        say(f"  {err['sample']}: {err['error']}")
say(f"Saved results: {OUT}")
say(
    "Next: review spot matching, QC distributions and annotation coverage before filtering or clustering."
)
report = "\n".join(summary)
(OUT / "Stage1B_Results.txt").write_text(report, encoding="utf-8")
print("\n" + report)

# ============================================================================
# Audit gene identifiers and establish the shared feature set
# ============================================================================

from pathlib import Path
import gzip
import io
import tarfile
import traceback
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = PROJECT_ROOT / "Data"
OUT = ROOT.parent / "Results" / "Stage1C"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
reports, stats, gene_sets, symbol_sets, errors, issues = [], [], {}, {}, [], []
spot_tables = []


def note(x=""):
    reports.append(str(x))


def member_stream(tar, path):
    raw = tar.extractfile(tar.getmember(path))
    if raw is None:
        raise ValueError(f"Cannot open {path}")
    buf = io.BufferedReader(raw)
    if buf.peek(2)[:2] == b"\x1f\x8b":
        return gzip.GzipFile(fileobj=buf)
    return buf


def read_tsv(tar, path):
    with member_stream(tar, path) as f:
        return pd.read_csv(f, sep="\t", header=None, dtype=str, keep_default_na=False)


def quartile(a):
    if len(a) == 0:
        return "NA"
    p = np.percentile(a, [25, 50, 75])
    return f"median {p[1]:,.1f}; IQR [{p[0]:,.1f}, {p[2]:,.1f}]"


note("PROJECT 2 | STAGE 1C — ESSENTIAL QC AND GENE HARMONISATION")
note("Diagnostics only. No filtering, normalisation, or data alteration.")
note(
    "Mitochondrial fractions require MT- gene symbols; absent MT- features are flagged."
)
note("")

for filename in ["filtered_count_matrices.tar.gz", "spatial.tar.gz", "metadata.tar.gz"]:
    if not (ROOT / filename).is_file():
        raise FileNotFoundError(ROOT / filename)

with (
    tarfile.open(ROOT / "filtered_count_matrices.tar.gz", "r:gz") as ct,
    tarfile.open(ROOT / "spatial.tar.gz", "r:gz") as st,
    tarfile.open(ROOT / "metadata.tar.gz", "r:gz") as mt,
):
    for sample in SAMPLES:
        note("=" * 65)
        note(f"SAMPLE: {sample}")
        try:
            prefix = f"filtered_count_matrices/{sample}_filtered_count_matrix"
            features = read_tsv(ct, f"{prefix}/features.tsv.gz")
            bars = read_tsv(ct, f"{prefix}/barcodes.tsv.gz").iloc[:, 0].astype(str)
            # Wu archive uses different feature-table formats. Some tables
            # contain ONLY gene symbols (one column), not Ensembl IDs.
            # Do not misrepresent symbols as stable gene identifiers.
            if features.shape[1] == 1:
                symbols = features.iloc[:, 0].astype(str)
                ids = symbols.copy()  # comparison key only; NOT an Ensembl ID
                feature_key_type = "gene_symbol_only"
            elif features.shape[1] >= 2:
                ids = features.iloc[:, 0].astype(str)
                symbols = features.iloc[:, 1].astype(str)
                feature_key_type = "feature_first_column"
            else:
                raise ValueError("Empty gene feature table")
            if ids.eq("").any() or symbols.eq("").any():
                raise ValueError("Empty gene comparison key or symbol")
            if ids.duplicated().any() or bars.duplicated().any():
                raise ValueError("Duplicate gene ID or barcode")
            with member_stream(ct, f"{prefix}/matrix.mtx.gz") as f:
                m = mmread(f)
            m = m.tocsr() if sparse.issparse(m) else sparse.csr_matrix(m)
            m.sum_duplicates()
            if m.shape != (len(features), len(bars)):
                raise ValueError(f"Incompatible matrix shape: {m.shape}")
            if (m.data < 0).any():
                raise ValueError("Negative counts")

            totals = np.asarray(m.sum(axis=0)).ravel()
            detected = np.asarray(m.getnnz(axis=0)).ravel()
            mt_mask = symbols.str.upper().str.startswith("MT-").to_numpy()
            mito = (
                np.asarray(m[mt_mask, :].sum(axis=0)).ravel()
                if mt_mask.any()
                else np.zeros(len(bars))
            )
            mito_pct = np.divide(
                mito * 100, totals, out=np.full(len(bars), np.nan), where=totals > 0
            )
            with member_stream(
                st, f"spatial/{sample}_spatial/tissue_positions_list.csv"
            ) as f:
                pos = pd.read_csv(
                    f,
                    header=None,
                    names=[
                        "barcode",
                        "in_tissue",
                        "array_row",
                        "array_col",
                        "pixel_row",
                        "pixel_col",
                    ],
                )
            with member_stream(mt, f"metadata/{sample}_metadata.csv") as f:
                meta = pd.read_csv(f, index_col=0)
            pos = pos.set_index("barcode")
            pos.index = pos.index.astype(str)
            meta.index = meta.index.astype(str)
            if pos.index.has_duplicates or meta.index.has_duplicates:
                raise ValueError("Duplicated position or metadata barcode")
            missing_pos = int((~bars.isin(pos.index)).sum())
            missing_meta = int((~bars.isin(meta.index)).sum())
            if missing_pos or missing_meta:
                raise ValueError(
                    f"Missing barcodes in positions/metadata: {missing_pos}/{missing_meta}"
                )
            pos = pos.reindex(bars)
            meta = meta.reindex(bars)
            classification = meta["Classification"].astype("string")
            out_tissue = (
                pd.to_numeric(pos["in_tissue"], errors="coerce").ne(1).to_numpy()
            )
            artifact = classification.str.contains(
                "artefact|artifact", case=False, na=False
            ).to_numpy(dtype=bool)
            missing_path = int(classification.isna().sum())
            zero_spots = int((totals == 0).sum())
            high_mito = int(np.sum(mito_pct > 20)) if mt_mask.any() else 0
            ids_set = set(ids)
            symbols_set = set(symbols[symbols.ne("")])
            gene_sets[sample] = ids_set
            symbol_sets[sample] = symbols_set
            frame = pd.DataFrame(
                {
                    "sample": sample,
                    "barcode": bars.to_numpy(),
                    "total_counts": totals,
                    "genes_detected": detected,
                    "pct_mito": mito_pct,
                    "in_tissue": pos["in_tissue"].to_numpy(),
                    "classification": classification.fillna("[MISSING]").to_numpy(),
                    "is_artefact": artifact,
                }
            )
            spot_tables.append(frame)
            abnormal = frame[
                (frame["in_tissue"] != 1)
                | frame["is_artefact"]
                | (frame["total_counts"] == 0)
            ]
            if not abnormal.empty:
                issues.append(abnormal)
            stats.append(
                {
                    "sample": sample,
                    "feature_key_type": feature_key_type,
                    "subtype": ";".join(
                        sorted(meta["subtype"].dropna().astype(str).unique())
                    ),
                    "genes": len(ids),
                    "spots": len(bars),
                    "median_counts": float(np.median(totals)),
                    "median_genes": float(np.median(detected)),
                    "mito_gene_features": int(mt_mask.sum()),
                    "median_pct_mito": float(np.nanmedian(mito_pct))
                    if np.isfinite(mito_pct).any()
                    else np.nan,
                    "pct_mito_gt20_n": high_mito,
                    "zero_count_spots": zero_spots,
                    "out_of_tissue_spots": int(out_tissue.sum()),
                    "artefact_spots": int(artifact.sum()),
                    "unannotated_spots": missing_path,
                }
            )
            note(
                f"Spots {len(bars):,}; genes {len(ids):,}; unique gene symbols {len(symbols_set):,}; key format: {feature_key_type}"
            )
            note(f"Counts/spot: {quartile(totals)}")
            note(f"Detected genes/spot: {quartile(detected)}")
            note(f"Mitochondrial genes detected in annotation: {int(mt_mask.sum())}")
            note(
                f"Mitochondrial %: {quartile(mito_pct[np.isfinite(mito_pct)]) if mt_mask.any() else 'NOT AVAILABLE'}; >20%: {high_mito}"
            )
            note(
                f"Zero-count spots: {zero_spots}; out-of-tissue: {int(out_tissue.sum())}; artefact: {int(artifact.sum())}; missing pathology: {missing_path}"
            )
            note("No spots excluded.")
        except Exception as exc:
            msg = f"{sample}: {type(exc).__name__}: {exc}"
            errors.append(msg)
            note("ERROR: " + msg)
            note(traceback.format_exc())

note("\n" + "=" * 65)
note("GENE HARMONISATION — FEATURE-KEY MATCHING (ID WHERE AVAILABLE, OTHERWISE SYMBOL)")
if gene_sets:
    shared = (
        set.intersection(*gene_sets.values())
        if len(gene_sets) == len(SAMPLES)
        else set()
    )
    union = set.union(*gene_sets.values())
    shared_symbols = (
        set.intersection(*symbol_sets.values())
        if len(symbol_sets) == len(SAMPLES)
        else set()
    )
    note(f"Unique feature keys (union): {len(union):,}")
    note(
        f"Feature keys in ALL {len(gene_sets)} successfully loaded specimens: {len(shared):,}"
    )
    note(f"Gene symbols common to ALL {len(symbol_sets)}: {len(shared_symbols):,}")
    note(
        "Feature key formats: "
        + "; ".join(f"{s['sample']}={s['feature_key_type']}" for s in stats)
    )
    note(
        "CAUTION: one-column feature files contain gene symbols, not stable Ensembl IDs."
    )
    note(
        "Cross-specimen overlap by gene symbol is provisional until feature key formats are examined."
    )
    for sample in gene_sets:
        note(
            f"  {sample}: IDs {len(gene_sets[sample]):,}; shared IDs {len(gene_sets[sample] & shared):,}; sample-only vs union of all others {len(gene_sets[sample] - set.union(*(s for other, s in gene_sets.items() if other != sample))):,}"
        )
    pd.DataFrame({"feature_key": sorted(shared)}).to_csv(
        OUT / "shared_feature_keys_all_samples.csv", index=False
    )
    pd.DataFrame({"feature_key": sorted(union)}).to_csv(
        OUT / "all_feature_keys_union.csv", index=False
    )
    pd.DataFrame({"gene_symbol": sorted(shared_symbols)}).to_csv(
        OUT / "shared_gene_symbols_all_samples.csv", index=False
    )
    pairwise = []
    for i, a in enumerate(gene_sets):
        for b in list(gene_sets)[i + 1 :]:
            inter = len(gene_sets[a] & gene_sets[b])
            pairwise.append(
                {
                    "sample_a": a,
                    "sample_b": b,
                    "shared_feature_keys": inter,
                    "jaccard": inter / len(gene_sets[a] | gene_sets[b]),
                }
            )
    pd.DataFrame(pairwise).to_csv(OUT / "pairwise_gene_overlap.csv", index=False)
else:
    note("No successful samples: gene overlap unavailable.")

if stats:
    stat_df = pd.DataFrame(stats)
    stat_df.to_csv(OUT / "sample_essential_qc.csv", index=False)
    note("\n" + "=" * 65)
    note("FINAL CONSOLIDATED QC TABLE")
    note(stat_df.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    note(f"Spots inspected: {int(stat_df['spots'].sum()):,}")
    note(f"Out-of-tissue spots: {int(stat_df['out_of_tissue_spots'].sum())}")
    note(f"Artefact-labeled spots: {int(stat_df['artefact_spots'].sum())}")
    note(f"Missing pathology: {int(stat_df['unannotated_spots'].sum())}")
    pd.concat(spot_tables, ignore_index=True).to_csv(
        OUT / "spot_essential_qc_unfiltered.csv", index=False
    )
    if issues:
        pd.concat(issues, ignore_index=True).to_csv(
            OUT / "spots_for_review.csv", index=False
        )
if errors:
    note("\nERRORS")
    for e in errors:
        note(e)
note(f"\nCompleted: {len(stats)}/{len(SAMPLES)} specimens; raw data unchanged")
note(f"Outputs: {OUT}")
note(
    "NEXT: review these diagnostics before deciding any exclusions and subsequent normalisation."
)
report = "\n".join(reports)
(OUT / "Stage1C_Results.txt").write_text(report, encoding="utf-8")
print("\n" + report)

# ============================================================================
# Apply exclusions and write harmonised sparse matrices
# ============================================================================

from pathlib import Path
import tarfile, gzip, io, traceback
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = PROJECT_ROOT / "Data"
OUT = ROOT.parent / "Results" / "Stage1D"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
report = []


def log(s=""):
    report.append(str(s))


def fileobj(tf, path):
    raw = tf.extractfile(tf.getmember(path))
    if raw is None:
        raise ValueError(f"Cannot read {path}")
    buf = io.BufferedReader(raw)
    return gzip.GzipFile(fileobj=buf) if buf.peek(2)[:2] == b"\x1f\x8b" else buf


def read_symbols(tf, sample):
    p = f"filtered_count_matrices/{sample}_filtered_count_matrix/features.tsv.gz"
    with fileobj(tf, p) as f:
        df = pd.read_csv(f, sep="\t", header=None, dtype=str, keep_default_na=False)
    sym = df.iloc[:, 0] if df.shape[1] == 1 else df.iloc[:, 1]
    sym = sym.astype(str)
    if sym.eq("").any() or sym.duplicated().any():
        raise ValueError(f"{sample}: missing or duplicate symbols")
    return sym


def initial_load(tf, sample):
    p = f"filtered_count_matrices/{sample}_filtered_count_matrix"
    sym = read_symbols(tf, sample)
    with fileobj(tf, f"{p}/barcodes.tsv.gz") as f:
        bars = pd.read_csv(f, sep="\t", header=None, dtype=str).iloc[:, 0].astype(str)
    with fileobj(tf, f"{p}/matrix.mtx.gz") as f:
        m = mmread(f)
    m = m.tocsr() if sparse.issparse(m) else sparse.csr_matrix(m)
    m.sum_duplicates()
    if m.shape != (len(sym), len(bars)):
        raise ValueError(f"{sample}: matrix dimensions invalid {m.shape}")
    if bars.duplicated().any():
        raise ValueError(f"{sample}: duplicate barcodes")
    return sym, bars, m


log("PROJECT 2 | STAGE 1D — ANALYSIS DATA PREPARATION")
log(
    "Rules: exclude pathology Artefact and out-of-tissue spots; retain high-mito and missing-pathology spots."
)
log("Expression: shared 18,714 symbols (if verified); sparse unnormalised RAW counts.")
log("")
for name in ["filtered_count_matrices.tar.gz", "spatial.tar.gz", "metadata.tar.gz"]:
    if not (ROOT / name).is_file():
        raise FileNotFoundError(ROOT / name)

stats = []
failures = []
try:
    with (
        tarfile.open(ROOT / "filtered_count_matrices.tar.gz", "r:gz") as counts,
        tarfile.open(ROOT / "spatial.tar.gz", "r:gz") as spat,
        tarfile.open(ROOT / "metadata.tar.gz", "r:gz") as met,
    ):
        # Shared-gene list is defined from all six specimens independently of QC filtering.
        all_sym = {s: read_symbols(counts, s) for s in SAMPLES}
        shared = set.intersection(*(set(x) for x in all_sym.values()))
        common = sorted(shared)
        if len(common) < 15000:
            raise ValueError(f"Unexpected shared gene count: {len(common)}")
        pd.DataFrame({"gene_symbol": common}).to_csv(
            OUT / "shared_gene_symbols_ordered.csv", index=False
        )
        log(f"Common gene symbols across all six: {len(common):,}")
        for sample in SAMPLES:
            try:
                sym, bars, m = initial_load(counts, sample)
                with fileobj(
                    spat, f"spatial/{sample}_spatial/tissue_positions_list.csv"
                ) as f:
                    pos = pd.read_csv(
                        f,
                        header=None,
                        names=[
                            "barcode",
                            "in_tissue",
                            "array_row",
                            "array_col",
                            "pixel_row",
                            "pixel_col",
                        ],
                    ).set_index("barcode")
                with fileobj(met, f"metadata/{sample}_metadata.csv") as f:
                    meta = pd.read_csv(f, index_col=0)
                pos.index = pos.index.astype(str)
                meta.index = meta.index.astype(str)
                if pos.index.has_duplicates or meta.index.has_duplicates:
                    raise ValueError("Duplicate spatial/metadata barcode")
                if not bars.isin(pos.index).all() or not bars.isin(meta.index).all():
                    raise ValueError("Barcode matching failed")
                pos = pos.reindex(bars)
                meta = meta.reindex(bars)
                original_totals = np.asarray(m.sum(axis=0)).ravel()
                pathology = meta["Classification"].astype("string")
                is_art = pathology.str.contains(
                    "artefact|artifact", case=False, na=False
                ).to_numpy(dtype=bool)
                in_tissue = (
                    pd.to_numeric(pos["in_tissue"], errors="coerce").eq(1).to_numpy()
                )
                keep = in_tissue & ~is_art & (original_totals > 0)
                excluded = int((~keep).sum())
                overlap = int(np.sum(~in_tissue & is_art))
                idx = pd.Index(sym)
                order = idx.get_indexer(common)
                if (order < 0).any():
                    raise ValueError("Missing common genes unexpectedly")
                # Matrix is genes x spots in archive; exported as spots x genes.
                X = m[order, :][:, keep].transpose().tocsr()
                if X.shape != (int(keep.sum()), len(common)):
                    raise ValueError("Output shape mismatch")
                if X.nnz and np.any(X.data < 0):
                    raise ValueError("Negative counts")
                kept_bars = bars[keep].to_numpy()
                obs = meta.iloc[np.flatnonzero(keep)].copy()
                xyz = pos.iloc[np.flatnonzero(keep)].copy()
                obs.index.name = "barcode"
                for c in xyz.columns:
                    obs[c] = xyz[c].to_numpy()
                obs["sample"] = sample
                obs["total_counts_full_panel"] = original_totals[keep]
                obs["total_counts_shared_genes"] = np.asarray(X.sum(axis=1)).ravel()
                obs["detected_genes_shared"] = np.asarray(X.getnnz(axis=1)).ravel()
                # Mito fractions must be calculated using the full original gene panel.
                mt_idx = sym.str.upper().str.startswith("MT-").to_numpy()
                mt_counts = np.asarray(m[mt_idx, :].sum(axis=0)).ravel()[keep]
                obs["pct_mito_full_panel"] = 100 * mt_counts / original_totals[keep]
                obs["Classification"] = (
                    pathology.iloc[np.flatnonzero(keep)].fillna("[MISSING]").to_numpy()
                )
                if not np.array_equal(obs.index.to_numpy(), kept_bars):
                    raise ValueError("Spot ordering mismatch")
                sparse.save_npz(
                    OUT / f"{sample}_shared_raw_spots_x_genes.npz", X, compressed=True
                )
                obs.to_csv(OUT / f"{sample}_spot_metadata.csv")
                stat = dict(
                    sample=sample,
                    initial_spots=m.shape[1],
                    removed_artifact=int(is_art.sum()),
                    removed_out_of_tissue=int((~in_tissue).sum()),
                    overlap_artifact_and_outside=overlap,
                    removed_union=excluded,
                    retained_spots=int(keep.sum()),
                    shared_genes=len(common),
                    missing_pathology_kept=int(
                        obs["Classification"].eq("[MISSING]").sum()
                    ),
                    retained_mito_gt20=int((obs["pct_mito_full_panel"] > 20).sum()),
                )
                stats.append(stat)
                log(
                    f"{sample}: {stat['initial_spots']:,} initial -> {stat['retained_spots']:,} retained; removed {excluded} (artefact={stat['removed_artifact']}, outside={stat['removed_out_of_tissue']}, overlap={overlap}); missing pathology kept={stat['missing_pathology_kept']}; mito>20% kept={stat['retained_mito_gt20']}"
                )
            except Exception as e:
                failures.append(f"{sample}: {type(e).__name__}: {e}")
                log("ERROR: " + failures[-1])
                log(traceback.format_exc())
except Exception as e:
    failures.append(f"GLOBAL: {type(e).__name__}: {e}")
    log("ERROR: " + failures[-1])
    log(traceback.format_exc())

log("\n" + "=" * 75 + "\nFINAL CONSOLIDATED RESULTS\n" + "=" * 75)
if stats:
    df = pd.DataFrame(stats)
    df.to_csv(OUT / "Stage1D_sample_summary.csv", index=False)
    log(df.to_string(index=False))
    log(f"Total initial spots: {df['initial_spots'].sum():,}")
    log(f"Total retained spots: {df['retained_spots'].sum():,}")
    log(f"Total unique exclusions: {df['removed_union'].sum():,}")
    log(f"Shared genes: {df['shared_genes'].iloc[0]:,}")
log(f"Completed: {len(stats)}/{len(SAMPLES)} samples")
log("All original archives remain unchanged; no normalisation applied.")
log(f"Outputs: {OUT}")
if failures:
    log("FAILURES:")
    for f in failures:
        log(f)
else:
    log("NEXT: review filtering and spot retention before normalisation / PCA.")
text = "\n".join(report)
(OUT / "Stage1D_Results.txt").write_text(text, encoding="utf-8")
print("\n" + text)
if failures:
    raise RuntimeError(f"Data preparation incomplete: {len(failures)} sample(s) failed; see Stage1D_Results.txt")
