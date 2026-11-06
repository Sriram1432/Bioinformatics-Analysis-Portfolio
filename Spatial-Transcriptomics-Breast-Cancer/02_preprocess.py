"""Normalise expression values and derive the retained PCA representations.

Calculations are performed independently within each specimen. The primary
representation follows shared-gene CP10K normalisation, log1p transformation,
variable-gene selection and deterministic PCA.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ============================================================================
# Within-specimen CP10K normalisation and log1p transformation
# ============================================================================

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import traceback


PROJECT = PROJECT_ROOT
SOURCE = PROJECT / "Results" / "Stage1D"
DEST = PROJECT / "Results" / "Stage2A"
DEST.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
TARGET = 10000.0
report = []
summary = []
failures = []


def log(x=""):
    report.append(str(x))


def qstr(x):
    q = np.quantile(x, [0, 0.25, 0.5, 0.75, 1])
    return f"min={q[0]:,.2f}; Q1={q[1]:,.2f}; median={q[2]:,.2f}; Q3={q[3]:,.2f}; max={q[4]:,.2f}"


log("PROJECT 2 | STAGE 2A — LIBRARY-SIZE NORMALISATION + LOG1P")
log(
    "Method: divide each spot by its total counts within the 18,714 shared genes; multiply by 10,000; apply ln(1+x)."
)
log(
    "Each specimen remains separate; raw counts preserved. No correction for specimen effects."
)
log("")

try:
    genes = pd.read_csv(SOURCE / "shared_gene_symbols_ordered.csv")[
        "gene_symbol"
    ].astype(str)
    if genes.duplicated().any() or genes.isna().any():
        raise ValueError("Invalid or duplicated shared gene symbols")
    log(f"Shared gene symbols: {len(genes):,}")
except Exception as e:
    raise RuntimeError(f"Missing/invalid Stage1D gene list: {e}") from e

for sample in SAMPLES:
    try:
        matrix_path = SOURCE / f"{sample}_shared_raw_spots_x_genes.npz"
        metadata_path = SOURCE / f"{sample}_spot_metadata.csv"
        X = sparse.load_npz(matrix_path).tocsr()
        obs = pd.read_csv(metadata_path, index_col=0)
        if X.shape != (len(obs), len(genes)):
            raise ValueError(
                f"Matrix {X.shape} does not match metadata rows {len(obs)} and genes {len(genes)}"
            )
        if X.nnz and (np.min(X.data) < 0 or not np.all(np.isfinite(X.data))):
            raise ValueError("Invalid raw counts")
        totals = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        if np.any(totals <= 0):
            raise ValueError(
                f"{np.sum(totals <= 0)} spots with zero shared-gene counts: investigate before normalisation"
            )
        metadata_totals = pd.to_numeric(
            obs["total_counts_shared_genes"], errors="raise"
        ).to_numpy()
        if not np.allclose(totals, metadata_totals, atol=0.01, rtol=0):
            raise ValueError("Shared gene count totals disagree with Stage1D metadata")
        # Use float64 arithmetic until AFTER log1p to avoid loss of
        # precision from float32 sparse scaling. Save as float32 afterward.
        scale = TARGET / totals
        Y = sparse.diags(scale, format="csr") @ X.astype(np.float64)
        Y = Y.tocsr()
        Y.data = np.log1p(Y.data)
        # Float32 is a storage format, not the calculation precision.
        Y = Y.astype(np.float32)
        Y.sort_indices()
        if Y.nnz != X.nnz or not np.all(np.isfinite(Y.data)):
            raise ValueError("Normalised sparse matrix validation failed")
        # Independently validate sampled spots using mathematical values;
        # compare aligned gene indices, not possibly unsorted CSR order.
        rng = np.random.default_rng(20261101)
        selected = rng.choice(X.shape[0], size=min(50, X.shape[0]), replace=False)
        max_abs_error = 0.0
        for i in selected:
            raw = X.getrow(int(i)).copy()
            raw.sort_indices()
            actual = Y.getrow(int(i))
            expected = np.log1p(raw.data.astype(np.float64) * scale[int(i)])
            if not np.array_equal(raw.indices, actual.indices):
                raise ValueError(f"Gene index mismatch for spot index {i}")
            error = float(
                np.max(np.abs(actual.data.astype(np.float64) - expected), initial=0.0)
            )
            max_abs_error = max(max_abs_error, error)
            if not np.allclose(actual.data, expected, rtol=1e-6, atol=1e-6):
                raise ValueError(
                    f"Numeric validation failed for spot index {i}; max abs error={error:.3g}"
                )
        log(
            f"{sample} | validation: {len(selected)} spots checked; max abs error={max_abs_error:.3g}"
        )
        out_path = DEST / f"{sample}_log1p_CP10K_spots_x_genes.npz"
        sparse.save_npz(out_path, Y, compressed=True)
        record = {
            "sample": sample,
            "spots": X.shape[0],
            "shared_genes": X.shape[1],
            "nonzero_entries": X.nnz,
            "median_raw_shared_counts": float(np.median(totals)),
            "min_raw_shared_counts": float(totals.min()),
            "max_raw_shared_counts": float(totals.max()),
            "median_raw_full_panel_counts": float(
                pd.to_numeric(obs["total_counts_full_panel"]).median()
            ),
            "zero_count_spots": int(np.sum(totals == 0)),
            "nonfinite_normalized_values": int(np.sum(~np.isfinite(Y.data))),
            "output_MB": round(out_path.stat().st_size / (1024**2), 2),
        }
        summary.append(record)
        log(
            f"{sample} | PASS | {X.shape[0]:,} spots × {X.shape[1]:,} genes | nnz={X.nnz:,}"
        )
        log(f"  Raw counts (shared genes): {qstr(totals)}")
        log(
            f"  Log-normalised nonzero value range: [{Y.data.min():.4f}, {Y.data.max():.4f}]"
        )
        log(f"  Saved: {out_path.name} ({record['output_MB']:.2f} MB)")
        del X, Y
    except Exception as e:
        failures.append(f"{sample}: {type(e).__name__}: {e}")
        log("ERROR: " + failures[-1])
        log(traceback.format_exc())

log("\n" + "=" * 78)
log("FINAL ANALYSIS SUMMARY")
log("=" * 78)
if summary:
    df = pd.DataFrame(summary)
    df.to_csv(DEST / "Stage2A_sample_summary.csv", index=False)
    log(df.to_string(index=False))
    log(f"Total normalised spots: {df.spots.sum():,}")
    log(f"Genes per specimen: {len(genes):,}")
log(f"Successfully processed: {len(summary)}/{len(SAMPLES)}")
log("Raw Stage1D matrices: unchanged.")
log(f"Outputs: {DEST}")
if failures:
    log("ERRORS:")
    for e in failures:
        log(e)
    log("Do not proceed to Stage 2B until these errors are resolved.")
else:
    log(
        "NEXT: inspect sample-level variability and select informative genes for PCA (Stage 2B)."
    )
full_report = "\n".join(report)
(DEST / "Stage2A_Results.txt").write_text(full_report, encoding="utf-8")
print("\n" + full_report)

# ============================================================================
# Variable-gene selection and primary PCA
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA


PROJECT = PROJECT_ROOT
IN_DIR = PROJECT / "Results" / "Stage2A"
META_DIR = PROJECT / "Results" / "Stage1D"
OUT_DIR = PROJECT / "Results" / "Stage2B"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
MIN_DETECTION_FRACTION = 0.05
N_HVG = 2000
N_PCS = 30
SEED = 20261101

messages, results, errors = [], [], []


def log(line=""):
    messages.append(str(line))


log("PROJECT 2 | STAGE 2B — PER-SPECIMEN VARIABLE GENES + PCA")
log("Input: Stage2A log1p(CP10K); no correction for specimen effects.")
log("Methods: within-sample detection >= 5%; top 2,000 by log-expression variance;")
log(
    "center and z-score selected genes, cap z-scores at +/-10, run 30-component randomized PCA."
)
log("PCA represents transcriptional variation among MULTICELLULAR tissue spots.")
log("")

genes = (
    pd.read_csv(META_DIR / "shared_gene_symbols_ordered.csv")["gene_symbol"]
    .astype(str)
    .to_numpy()
)
if len(set(genes)) != len(genes):
    raise ValueError("Gene symbols are not unique")

for sample in SAMPLES:
    try:
        X = sparse.load_npz(IN_DIR / f"{sample}_log1p_CP10K_spots_x_genes.npz").tocsr()
        obs = pd.read_csv(META_DIR / f"{sample}_spot_metadata.csv", index_col=0)
        if X.shape != (len(obs), len(genes)) or X.shape[0] < 50:
            raise ValueError(
                f"Inconsistent dimensions: matrix={X.shape}, metadata={obs.shape}, genes={len(genes)}"
            )
        if X.data.size and not np.isfinite(X.data).all():
            raise ValueError("Nonfinite normalised expression values")
        n = X.shape[0]
        detect = np.asarray((X > 0).sum(axis=0)).ravel() / n
        means = np.asarray(X.mean(axis=0)).ravel().astype(np.float64)
        squared = X.copy().astype(np.float64)
        squared.data **= 2
        variance = np.maximum(0, np.asarray(squared.mean(axis=0)).ravel() - means**2)
        eligible = np.flatnonzero(
            (detect >= MIN_DETECTION_FRACTION) & (variance > 1e-10)
        )
        if len(eligible) < N_HVG:
            raise ValueError(
                f"Only {len(eligible)} genes meet detection and variance requirements"
            )
        # Stable ranking under ties; choose highest variance, then alphabetical gene name.
        rank_order = np.lexsort((genes[eligible], -variance[eligible]))
        chosen = eligible[rank_order[:N_HVG]]
        chosen = chosen[np.argsort(genes[chosen])]
        hvg = pd.DataFrame(
            {
                "gene_symbol": genes[chosen],
                "mean_log_expression": means[chosen],
                "variance_log_expression": variance[chosen],
                "fraction_spots_detected": detect[chosen],
            }
        )
        hvg.to_csv(OUT_DIR / f"{sample}_top{N_HVG}_variable_genes.csv", index=False)
        # 4,847 x 2,000 float32 <= 39 MB per sample. Per-sample centering avoids
        # conflating within-tumour spatial patterns with between-patient effects.
        Z = X[:, chosen].toarray().astype(np.float32)
        mu = Z.mean(axis=0, dtype=np.float64)
        sd = Z.std(axis=0, dtype=np.float64)
        if np.any(sd <= 0):
            raise ValueError("Zero-variance gene after selection")
        Z = np.clip((Z.astype(np.float64) - mu) / sd, -10, 10).astype(np.float32)
        pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=SEED)
        scores = pca.fit_transform(Z)
        if scores.shape != (n, N_PCS) or not np.isfinite(scores).all():
            raise ValueError("PCA failed numerical validation")
        pc_cols = [f"PC{i + 1}" for i in range(N_PCS)]
        coords = pd.DataFrame(scores, columns=pc_cols, index=obs.index)
        coords.index.name = "barcode"
        coords.to_csv(OUT_DIR / f"{sample}_PCA30_spot_scores.csv")
        ev = pd.DataFrame(
            {
                "component": pc_cols,
                "variance_ratio": pca.explained_variance_ratio_,
                "cumulative_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
            }
        )
        ev.to_csv(OUT_DIR / f"{sample}_PCA_explained_variance.csv", index=False)
        pd.DataFrame(
            pca.components_.T, index=genes[chosen], columns=pc_cols
        ).rename_axis("gene_symbol").to_csv(OUT_DIR / f"{sample}_PCA_gene_loadings.csv")
        rec = {
            "sample": sample,
            "spots": n,
            "genes_available": X.shape[1],
            "genes_detected_5pct": len(eligible),
            "HVG_selected": len(chosen),
            "PC1_variance_pct": round(100 * ev.variance_ratio.iloc[0], 2),
            "PC1to10_variance_pct": round(100 * ev.variance_ratio.iloc[:10].sum(), 2),
            "PC1to30_variance_pct": round(100 * ev.variance_ratio.sum(), 2),
        }
        results.append(rec)
        log(
            f"{sample}: PASS | {n:,} spots | eligible genes {len(eligible):,} | HVGs {len(chosen):,} | "
            f"PC1 {rec['PC1_variance_pct']:.2f}% | PC1–10 {rec['PC1to10_variance_pct']:.2f}% | "
            f"PC1–30 {rec['PC1to30_variance_pct']:.2f}%"
        )
        del X, Z, squared, scores
    except Exception as exc:
        errors.append(f"{sample}: {type(exc).__name__}: {exc}")
        log("ERROR: " + errors[-1])
        log(traceback.format_exc())

log("\n" + "=" * 78)
log("FINAL ANALYSIS SUMMARY")
log("=" * 78)
if results:
    table = pd.DataFrame(results)
    table.to_csv(OUT_DIR / "Stage2B_sample_summary.csv", index=False)
    log(table.to_string(index=False))
    log(f"Samples completed: {len(results)}/{len(SAMPLES)}")
    log(f"Total PCA spots: {sum(r['spots'] for r in results):,}")
    # Compare the selected variable gene sets to distinguish general from specimen-specific variability.
    gene_sets = [
        set(pd.read_csv(OUT_DIR / f"{s}_top{N_HVG}_variable_genes.csv")["gene_symbol"])
        for s in SAMPLES
        if s in {r["sample"] for r in results}
    ]
    if len(gene_sets) == len(SAMPLES):
        shared_hvg = set.intersection(*gene_sets)
        log(
            f"Variable genes shared across all six HVG sets: {len(shared_hvg)} / {N_HVG}"
        )
        pd.DataFrame({"gene_symbol": sorted(shared_hvg)}).to_csv(
            OUT_DIR / "HVG_shared_all_six.csv", index=False
        )
log("Raw matrices and Stage2A normalised matrices: unchanged.")
log("No UMAP, clusters, or cell identities have been inferred yet.")
log(f"Output directory: {OUT_DIR}")
if errors:
    log("ERRORS:")
    for err in errors:
        log(err)
    log("Resolve these errors before Stage 2C.")
else:
    log(
        "NEXT: evaluate PCA structure against tissue pathology, then consider UMAP/clustering."
    )
report = "\n".join(messages)
(OUT_DIR / "Stage2B_Results.txt").write_text(report, encoding="utf-8")
print("\n" + report)

# ============================================================================
# Descriptive PCA, pathology and spatial diagnostics
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = PROJECT_ROOT
META = ROOT / "Results" / "Stage1D"
PCA = ROOT / "Results" / "Stage2B"
OUT = ROOT / "Results" / "Stage2C"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
MIN_CATEGORY_N = 20  # for readable legends and descriptive category association; no filtering of spots
messages, rows, errors = [], [], []


def log(x=""):
    messages.append(str(x))


def eta_squared(x, group):
    # Descriptive fraction of between-annotation variance; NOT an independence-based significance test.
    df = pd.DataFrame(
        {"x": np.asarray(x, dtype=float), "g": pd.Series(group).to_numpy()}
    ).dropna()
    if df.empty:
        return np.nan
    grand = df.x.mean()
    ss_all = np.square(df.x - grand).sum()
    if ss_all <= 0:
        return np.nan
    ss_between = sum(
        len(a) * float((a.x.mean() - grand) ** 2) for _, a in df.groupby("g")
    )
    return ss_between / ss_all


def save_plot(fig, dest):
    fig.savefig(dest, dpi=230, bbox_inches="tight", format="webp")
    plt.close(fig)


log("PROJECT 2 | STAGE 2C — PCA–PATHOLOGY AND SPATIAL DIAGNOSTICS")
log(
    "Descriptive only; pathology and spatial spots are NOT statistically independent replicates."
)
log("No clustering, cell-type assignment, filtering or new normalisation.")
log("")
for sample in SAMPLES:
    try:
        obs = pd.read_csv(META / f"{sample}_spot_metadata.csv", index_col="barcode")
        pc = pd.read_csv(PCA / f"{sample}_PCA30_spot_scores.csv", index_col="barcode")
        ev = pd.read_csv(PCA / f"{sample}_PCA_explained_variance.csv")
        ld = pd.read_csv(
            PCA / f"{sample}_PCA_gene_loadings.csv", index_col="gene_symbol"
        )
        if not pc.index.equals(obs.index):
            raise ValueError(
                "Spot metadata and PCA row order differ; refusing positional alignment"
            )
        needed = [
            "Classification",
            "pixel_row",
            "pixel_col",
            "total_counts_shared_genes",
            "pct_mito_full_panel",
        ]
        missing = [c for c in needed if c not in obs]
        if missing:
            raise ValueError(f"Missing metadata columns: {missing}")
        labels = obs["Classification"].fillna("[MISSING]").astype(str)
        counts = labels.value_counts()
        labels_plot = labels.where(labels.map(counts) >= MIN_CATEGORY_N, "Other (n<20)")
        etas = {f"PC{i}": eta_squared(pc[f"PC{i}"], labels_plot) for i in range(1, 11)}
        cor = {}
        for cov in ["total_counts_shared_genes", "pct_mito_full_panel"]:
            x = np.asarray(obs[cov], float)
            x = np.log1p(x) if cov == "total_counts_shared_genes" else x
            for axis in ["PC1", "PC2"]:
                cor[f"{axis}_vs_{cov}"] = float(
                    np.corrcoef(pc[axis].to_numpy(float), x)[0, 1]
                )
        # PC1 vs PC2 colored by annotations; no category assigned to a gene-expression cluster.
        categories = labels_plot.value_counts().index.tolist()
        fig, ax = plt.subplots(figsize=(8.5, 6.3))
        cmap = plt.get_cmap("tab20")
        for i, cat in enumerate(categories):
            mask = (labels_plot == cat).to_numpy()
            ax.scatter(
                pc.loc[mask, "PC1"],
                pc.loc[mask, "PC2"],
                s=7,
                alpha=0.58,
                color=cmap(i % 20),
                label=f"{cat} (n={mask.sum()})",
                rasterized=True,
            )
        ax.set_xlabel(f"PC1 ({100 * ev.variance_ratio.iloc[0]:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({100 * ev.variance_ratio.iloc[1]:.1f}% variance)")
        ax.set_title(f"{sample} | PCA annotated by pathology")
        ax.legend(
            title="Histopathology",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=8,
            markerscale=1.7,
        )
        fig.tight_layout()
        save_plot(fig, OUT / f"{sample}_PCA_pathology.webp")
        # Spatial projection of PC1; image coordinates (pixel_col=x, pixel_row=y).
        fig, ax = plt.subplots(figsize=(7.2, 6.0))
        scat = ax.scatter(
            obs.pixel_col,
            obs.pixel_row,
            c=pc.PC1,
            s=7,
            cmap="coolwarm",
            alpha=0.85,
            rasterized=True,
        )
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Full-resolution image x (pixels)")
        ax.set_ylabel("Full-resolution image y (pixels)")
        ax.set_title(f"{sample} | spatial distribution of PC1")
        fig.colorbar(scat, ax=ax, label="PC1 score")
        fig.tight_layout()
        save_plot(fig, OUT / f"{sample}_Spatial_PC1.webp")
        # Top absolute loadings (sign arbitrary, values describe correlation-like PCA directions).
        loadings = ld["PC1"].sort_values()
        neg = loadings.head(10)
        pos = loadings.tail(10).sort_values(ascending=False)
        rec = {
            "sample": sample,
            "spots": len(obs),
            "pathology_categories": int(labels.nunique()),
            "largest_category": counts.index[0],
            "largest_category_n": int(counts.iloc[0]),
            "PC1_variance_pct": round(float(100 * ev.variance_ratio.iloc[0]), 2),
            "PC2_variance_pct": round(float(100 * ev.variance_ratio.iloc[1]), 2),
            "pathology_PC1_eta2_descriptive": round(float(etas["PC1"]), 4),
            "pathology_PC2_eta2_descriptive": round(float(etas["PC2"]), 4),
            **{k: round(v, 4) for k, v in cor.items()},
        }
        rows.append(rec)
        log(
            f"{sample}: PASS | spots={len(obs):,} | pathology groups={len(counts)} | "
            f"PC1={rec['PC1_variance_pct']:.2f}% | pathology eta² on PC1={rec['pathology_PC1_eta2_descriptive']:.3f}"
        )
        log(
            "  PC1 positive loadings: "
            + ", ".join(f"{k} ({v:+.3f})" for k, v in pos.items())
        )
        log(
            "  PC1 negative loadings: "
            + ", ".join(f"{k} ({v:+.3f})" for k, v in neg.items())
        )
        log(
            "  Pathology category counts: "
            + "; ".join(f"{k}={v}" for k, v in counts.items())
        )
        log(
            f"  Correlation PC1 with log library size: {cor['PC1_vs_total_counts_shared_genes']:+.3f}; "
            f"with mitochondrial fraction: {cor['PC1_vs_pct_mito_full_panel']:+.3f}"
        )
        pd.DataFrame({"category": counts.index, "n_spots": counts.values}).to_csv(
            OUT / f"{sample}_pathology_counts.csv", index=False
        )
        pd.DataFrame(
            {"PC": list(etas), "eta2_descriptive": list(etas.values())}
        ).to_csv(OUT / f"{sample}_PC_pathology_eta2.csv", index=False)
    except Exception as e:
        errors.append(f"{sample}: {type(e).__name__}: {e}")
        log("ERROR: " + errors[-1])
        log(traceback.format_exc())

log("\n" + "=" * 82)
log("FINAL ANALYSIS SUMMARY")
log("=" * 82)
if rows:
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "Stage2C_PCA_pathology_summary.csv", index=False)
    log(frame.to_string(index=False))
log(f"Samples completed: {len(rows)}/{len(SAMPLES)}")
log(f"Saved WebP figures: {2 * len(rows)}")
log(
    "NOTE: eta² is descriptive; spot-level significance tests would be pseudoreplicated."
)
log(
    "NOTE: PC axes/signs are sample-specific; do not compare their signs across tumours."
)
log(
    "NOTE: Pathological composition, library depth and technical effects can all contribute to PC variation."
)
if errors:
    log("ERRORS:")
    for error in errors:
        log(error)
else:
    log(
        "NEXT: review PCA pathology structure and potential QC confounding before UMAP/clustering."
    )
log(f"Outputs: {OUT}")
report = "\n".join(messages)
(OUT / "Stage2C_Results.txt").write_text(report, encoding="utf-8")
print(report)

# ============================================================================
# Sequencing-depth sensitivity comparison
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA


ROOT = PROJECT_ROOT
RAW = ROOT / "Results" / "Stage1D"
OLD = ROOT / "Results" / "Stage2B"
OUT = ROOT / "Results" / "Stage2D"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
PCS = 15
SEED = 20261101
MIN_GROUP = 20
report, summary, errors = [], [], []


def log(s=""):
    report.append(str(s))


def eta2(values, groups):
    v = np.asarray(values, dtype=np.float64)
    g = np.asarray(groups, dtype=str)
    good = np.isfinite(v) & (g != "[MISSING]")
    v, g = v[good], g[good]
    if len(v) < 3 or len(np.unique(g)) < 2:
        return np.nan
    mu = v.mean()
    tot = np.square(v - mu).sum()
    if tot <= 0:
        return np.nan
    return float(
        sum(
            int(np.sum(g == label)) * (v[g == label].mean() - mu) ** 2
            for label in np.unique(g)
        )
        / tot
    )


def correlation(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


log("PROJECT 2 | STAGE 2D — DEPTH ASSOCIATION: PCA SENSITIVITY CHECK")
log("Existing: log1p(CP10K), gene z-score + clipping, PCA (Stage 2B).")
log("Alternative: sqrt(CP10K) from RAW counts, gene-centering only, PCA.")
log("Same 2,000 HVGs per specimen in both comparisons; labels only for diagnostics.")
log("NO automatic regression, filtering, inferential tests or clustering.")

all_genes = (
    pd.read_csv(RAW / "shared_gene_symbols_ordered.csv")["gene_symbol"]
    .astype(str)
    .to_numpy()
)
if len(all_genes) != len(set(all_genes)):
    raise ValueError("Duplicate shared gene names")
lookup = pd.Index(all_genes)

for sample in SAMPLES:
    try:
        obs = pd.read_csv(RAW / f"{sample}_spot_metadata.csv", index_col="barcode")
        old = pd.read_csv(OLD / f"{sample}_PCA30_spot_scores.csv", index_col="barcode")
        hvg = (
            pd.read_csv(OLD / f"{sample}_top2000_variable_genes.csv")["gene_symbol"]
            .astype(str)
            .to_numpy()
        )
        X = sparse.load_npz(RAW / f"{sample}_shared_raw_spots_x_genes.npz").tocsr()
        if X.shape != (len(obs), len(all_genes)) or not obs.index.equals(old.index):
            raise ValueError(
                "Matrix dimensions or barcode order differ from previous stages"
            )
        col_idx = lookup.get_indexer(hvg)
        if np.any(col_idx < 0) or len(np.unique(col_idx)) != 2000:
            raise ValueError("Missing or duplicate HVG names")
        totals = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        if np.any(totals <= 0):
            raise ValueError("Zero raw-count spots found")
        depth = np.log1p(totals)
        # Work only with the 2,000 HVGs; peak dense array ~80 MB.
        Y = X[:, col_idx].toarray().astype(np.float64)
        Y *= (1e4 / totals)[:, None]
        np.sqrt(Y, out=Y)
        Y -= Y.mean(axis=0)
        # Do not standardise each gene: avoid inflating low-abundance genes.
        model = PCA(n_components=PCS, svd_solver="randomized", random_state=SEED)
        alt = model.fit_transform(Y)
        if not np.isfinite(alt).all():
            raise ValueError("Nonfinite alternative PCA scores")
        cols = [f"PC{i}" for i in range(1, PCS + 1)]
        scores = pd.DataFrame(alt, index=obs.index, columns=cols)
        scores.index.name = "barcode"
        scores.to_csv(OUT / f"{sample}_sqrtCP10K_PCA15_scores.csv")
        pd.DataFrame(
            {
                "gene_symbol": hvg,
                **{f"{p}_loading": model.components_[i, :] for i, p in enumerate(cols)},
            }
        ).to_csv(OUT / f"{sample}_sqrtCP10K_PCA15_loadings.csv", index=False)
        labels = obs["Classification"].fillna("[MISSING]").astype(str)
        counts = labels.value_counts()
        valid_groups = labels.where(labels.map(counts) >= MIN_GROUP, "Other_n_lt_20")
        valid_groups = valid_groups.where(labels != "[MISSING]", "[MISSING]")
        # Both PCA directions have arbitrary sign; use absolute r to compare.
        # PC1 may represent different biological directions across methods.
        old_depth = correlation(old.PC1, depth)
        alt_depth = correlation(scores.PC1, depth)
        old_path = eta2(old.PC1, valid_groups)
        alt_path = eta2(scores.PC1, valid_groups)
        rec = dict(
            sample=sample,
            spots=len(obs),
            old_PC1_abs_r_depth=round(abs(old_depth), 4),
            sqrt_PC1_abs_r_depth=round(abs(alt_depth), 4),
            old_PC1_pathology_eta2=round(old_path, 4),
            sqrt_PC1_pathology_eta2=round(alt_path, 4),
            old_PC2_abs_r_depth=round(abs(correlation(old.PC2, depth)), 4),
            sqrt_PC2_abs_r_depth=round(abs(correlation(scores.PC2, depth)), 4),
            sqrt_PC1_variance_pct=round(
                float(100 * model.explained_variance_ratio_[0]), 2
            ),
            sqrt_PC1to15_variance_pct=round(
                float(100 * model.explained_variance_ratio_.sum()), 2
            ),
        )
        summary.append(rec)
        # Each method's PC1 compared with depth AND all pathology groups.
        log(
            f"{sample}: PASS | old vs sqrt abs(r_depth)={rec['old_PC1_abs_r_depth']:.3f} / {rec['sqrt_PC1_abs_r_depth']:.3f} | "
            f"pathology eta²={old_path:.3f} / {alt_path:.3f}"
        )
        del X, Y, alt
    except Exception as exc:
        errors.append(f"{sample}: {type(exc).__name__}: {exc}")
        log("ERROR: " + errors[-1])
        log(traceback.format_exc())

log("\n" + "=" * 90)
log("FINAL ANALYSIS SUMMARY")
log("=" * 90)
if summary:
    frame = pd.DataFrame(summary)
    frame.to_csv(OUT / "Stage2D_PCA_method_comparison.csv", index=False)
    log(frame.to_string(index=False))
    log(f"Completed: {len(summary)}/{len(SAMPLES)} specimens")
    log(
        f"Median PC1 |r(depth)| existing: {frame.old_PC1_abs_r_depth.median():.3f}; sqrt: {frame.sqrt_PC1_abs_r_depth.median():.3f}"
    )
    log(
        f"Median PC1 pathology eta² existing: {frame.old_PC1_pathology_eta2.median():.3f}; sqrt: {frame.sqrt_PC1_pathology_eta2.median():.3f}"
    )
log(
    "Interpretation caution: pathology and depth can be biologically associated; lower |r| is not proof of superiority."
)
log(
    "PC1 axes can change biological meaning between methods; eta² is descriptive, not a significance test."
)
log(
    "Shared HVGs were originally selected by log-expression variance, not count-model residual dispersion."
)
log("No change to Stage2B or Stage2C outputs. No method selected automatically.")
if errors:
    log("ERRORS:")
    for item in errors:
        log(item)
else:
    log(
        "NEXT: inspect this comparison, then select PCA representation before neighbourhood graph/UMAP."
    )
log(f"Outputs: {OUT}")
result = "\n".join(report)
(OUT / "Stage2D_Results.txt").write_text(result, encoding="utf-8")
print(result)
if errors:
    raise RuntimeError(f"Preprocessing incomplete: {len(errors)} error(s); see Stage2D_Results.txt")
