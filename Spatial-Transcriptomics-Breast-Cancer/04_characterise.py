"""Characterise transcriptional clusters and curated expression programmes.

Marker rankings are descriptive within-specimen contrasts rather than formal
differential-expression tests. Programme scores are interpreted as molecular
signals in multicellular spots, not as direct cell-type assignments.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ============================================================================
# Descriptive cluster-associated gene ranking
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = PROJECT_ROOT
S1 = ROOT / "Results" / "Stage1D"
S2 = ROOT / "Results" / "Stage2A"
S3 = ROOT / "Results" / "Stage2E"
OUT = ROOT / "Results" / "Stage3A"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
TOP_N = 12
MIN_DETECTION_CLUSTER = 0.15
MIN_DETECTION_DIFFERENCE = 0.05
MIN_MEAN_LOG_DIFFERENCE = 0.25
report = []
results = []
top_rows = []
overview = []
errors = []


def log(msg=""):
    report.append(str(msg))


def analyze(sample, symbols):
    X = sparse.load_npz(S2 / f"{sample}_log1p_CP10K_spots_x_genes.npz").tocsr()
    obs = pd.read_csv(S1 / f"{sample}_spot_metadata.csv", index_col="barcode")
    assignment = pd.read_csv(S3 / f"{sample}_log_zscore_spot_clusters_UMAP.csv")
    if "barcode" not in assignment or "cluster" not in assignment:
        raise ValueError("Stage2E cluster assignment missing barcode/cluster columns")
    if assignment["barcode"].duplicated().any() or obs.index.has_duplicates:
        raise ValueError("Duplicate spot barcodes detected")
    if not np.array_equal(
        assignment["barcode"].astype(str).to_numpy(), obs.index.astype(str).to_numpy()
    ):
        raise ValueError(
            "Cluster and metadata spot order mismatch; refusing to infer order"
        )
    if X.shape != (len(obs), len(symbols)):
        raise ValueError(
            f"Expression shape {X.shape} != expected {(len(obs), len(symbols))}"
        )
    if X.nnz and (not np.isfinite(X.data).all() or (X.data < 0).any()):
        raise ValueError("Invalid normalised matrix values")
    labels = pd.to_numeric(assignment["cluster"], errors="raise").to_numpy(
        dtype=np.int64
    )
    groups, inverse = np.unique(labels, return_inverse=True)
    if len(groups) < 2:
        raise ValueError("Need at least two clusters")
    ng = len(groups)
    ns = X.shape[0]
    # Sparse membership multiplication calculates per-cluster gene sums without densifying X.
    membership = sparse.csr_matrix(
        (np.ones(ns, dtype=np.float64), (inverse, np.arange(ns))), shape=(ng, ns)
    )
    X64 = X.astype(np.float64)
    sums = (membership @ X64).toarray()
    detected = X.copy()
    detected.data = np.ones_like(detected.data, dtype=np.float32)
    detections = (membership @ detected).toarray()
    total_sums = sums.sum(axis=0)
    total_detected = detections.sum(axis=0)
    sizes = np.bincount(inverse, minlength=ng)
    mean = sums / sizes[:, None]
    fraction = detections / sizes[:, None]
    rest_mean = (total_sums[None, :] - sums) / (ns - sizes)[:, None]
    rest_fraction = (total_detected[None, :] - detections) / (ns - sizes)[:, None]
    delta = mean - rest_mean
    delta_detected = fraction - rest_fraction
    pathway_hint = {
        "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2", "CD24", "MUC1"],
        "immune": ["PTPRC", "CD3D", "CD3E", "CD79A", "MS4A1", "IGKC", "CD74"],
        "stromal_ECM": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "ACTA2"],
        "vascular": ["PECAM1", "VWF", "KDR", "EMCN", "CLDN5"],
        "cycling": ["MKI67", "TOP2A", "STMN1", "PCNA", "UBE2C"],
    }
    gene_to_idx = {gene: i for i, gene in enumerate(symbols)}
    guide_rows = []
    for kind, markers in pathway_hint.items():
        found = [gene_to_idx[g] for g in markers if g in gene_to_idx]
        for k, cl in enumerate(groups):
            guide_rows.append(
                dict(
                    sample=sample,
                    cluster=int(cl),
                    program=kind,
                    available_genes=len(found),
                    mean_log_expression=float(mean[k, found].mean())
                    if found
                    else np.nan,
                )
            )
    pd.DataFrame(guide_rows).to_csv(
        OUT / f"{sample}_candidate_signature_descriptive.csv", index=False
    )
    pathway = obs["Classification"].fillna("[MISSING]").astype(str).to_numpy()
    breakdown = pd.crosstab(
        pd.Series(labels, name="cluster"), pd.Series(pathway, name="pathology")
    )
    breakdown.to_csv(OUT / f"{sample}_cluster_pathology_counts.csv")
    marker_rows = []
    top_by_cluster = []
    for k, cl in enumerate(groups):
        eligibility = (
            (fraction[k] >= MIN_DETECTION_CLUSTER)
            & (delta_detected[k] >= MIN_DETECTION_DIFFERENCE)
            & (delta[k] >= MIN_MEAN_LOG_DIFFERENCE)
        )
        # Exclude mitochondrial and ribosomal genes from the primary marker shortlist,
        # but retain them in the complete table (all effects saved).
        priority = eligibility & np.array(
            [not str(g).startswith(("MT-", "RPS", "RPL")) for g in symbols]
        )
        ranked = np.flatnonzero(priority)
        ranked = ranked[np.argsort(-delta[k, ranked], kind="stable")]
        top_idx = ranked[:TOP_N]
        top_by_cluster.extend(list(top_idx[:2]))
        for rank, j in enumerate(top_idx, 1):
            row = dict(
                sample=sample,
                cluster=int(cl),
                rank=rank,
                gene_symbol=symbols[j],
                cluster_spots=int(sizes[k]),
                mean_log_cluster=round(float(mean[k, j]), 4),
                mean_log_rest=round(float(rest_mean[k, j]), 4),
                delta_mean_log=round(float(delta[k, j]), 4),
                detected_cluster_pct=round(float(100 * fraction[k, j]), 2),
                detected_rest_pct=round(float(100 * rest_fraction[k, j]), 2),
                delta_detected_pct=round(float(100 * delta_detected[k, j]), 2),
            )
            marker_rows.append(row)
            top_rows.append(row)
        overview.append(
            dict(
                sample=sample,
                cluster=int(cl),
                spots=int(sizes[k]),
                candidate_markers_passing=int(priority.sum()),
                top_markers=", ".join(symbols[j] for j in top_idx[:6]),
            )
        )
    pd.DataFrame(marker_rows).to_csv(
        OUT / f"{sample}_top_descriptive_markers.csv", index=False
    )
    # Save full gene-wise descriptive contrasts for transparency / downstream investigation.
    all_rows = []
    for k, cl in enumerate(groups):
        all_rows.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "cluster": int(cl),
                    "gene_symbol": symbols,
                    "mean_log_cluster": mean[k].astype(np.float32),
                    "mean_log_rest": rest_mean[k].astype(np.float32),
                    "delta_mean_log": delta[k].astype(np.float32),
                    "detected_cluster_frac": fraction[k].astype(np.float32),
                    "detected_rest_frac": rest_fraction[k].astype(np.float32),
                }
            )
        )
    pd.concat(all_rows, ignore_index=True).to_csv(
        OUT / f"{sample}_all_gene_descriptive_contrasts.csv.gz",
        index=False,
        compression="gzip",
    )
    # One standalone heatmap per specimen; group means are DESCRIPTIVE.
    uniq = []
    for j in top_by_cluster:
        if j not in uniq:
            uniq.append(j)
    if uniq:
        shown = mean[:, uniq]
        fig, ax = plt.subplots(
            figsize=(max(8, len(uniq) * 0.42), max(4.3, ng * 0.48 + 1.5))
        )
        im = ax.imshow(shown, aspect="auto", interpolation="nearest", cmap="viridis")
        ax.set_xticks(
            np.arange(len(uniq)),
            [symbols[j] for j in uniq],
            rotation=65,
            ha="right",
            fontsize=9,
        )
        ax.set_yticks(np.arange(ng), [str(int(x)) for x in groups])
        ax.set_xlabel("Cluster-associated gene", fontweight="bold")
        ax.set_ylabel("Transcriptional spot cluster", fontweight="bold")
        ax.set_title(f"{sample} | Mean log1p(CP10K) expression", fontweight="bold")
        fig.colorbar(im, ax=ax, label="Mean log-normalised expression", shrink=0.85)
        fig.tight_layout()
        fig.savefig(
            OUT / f"{sample}_cluster_marker_heatmap.webp", dpi=230, bbox_inches="tight"
        )
        plt.close(fig)
    dominant = breakdown.idxmax(axis=1)
    log(
        f"{sample}: PASS | {ns:,} spots | {ng} clusters | {len(marker_rows)} shortlisted markers"
    )
    current = [r for r in overview if r["sample"] == sample]
    for r in current:
        log(
            f"  Cluster {r['cluster']} (n={r['spots']}): {r['top_markers'] or 'no markers passing thresholds'}; "
            f"candidate markers={r['candidate_markers_passing']}"
        )
    results.append(
        dict(
            sample=sample,
            spots=ns,
            clusters=ng,
            shortlisted_marker_rows=len(marker_rows),
            marker_gene_union=len(set(x["gene_symbol"] for x in marker_rows)),
        )
    )


log("PROJECT 2 | STAGE 3A — DESCRIPTIVE SPATIAL-CLUSTER MARKER PROGRAMMES")
log("Primary: Stage2E log_zscore clusters; expression values: Stage2A log1p(CP10K).")
log(
    f"Marker shortlist: cluster detection >= {MIN_DETECTION_CLUSTER:.0%}, detection difference >= {MIN_DETECTION_DIFFERENCE:.0%}, mean-log difference >= {MIN_MEAN_LOG_DIFFERENCE:.2f}."
)
log(
    "Effects are descriptive (NOT differential-expression p-values); Visium spots contain multiple cells."
)
try:
    symbols = (
        pd.read_csv(S1 / "shared_gene_symbols_ordered.csv")["gene_symbol"]
        .astype(str)
        .tolist()
    )
    if len(symbols) != len(set(symbols)):
        raise ValueError("Duplicate shared gene symbols")
    for sample in SAMPLES:
        try:
            analyze(sample, symbols)
        except Exception as exc:
            errors.append(f"{sample}: {type(exc).__name__}: {exc}")
            log("ERROR: " + errors[-1])
            log(traceback.format_exc())
except Exception as exc:
    errors.append(f"GLOBAL: {type(exc).__name__}: {exc}")
    log("ERROR: " + errors[-1])
    log(traceback.format_exc())
log("\n" + "=" * 90 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 90)
if results:
    df = pd.DataFrame(results)
    df.to_csv(OUT / "Stage3A_sample_summary.csv", index=False)
    log(df.to_string(index=False))
if overview:
    ov = pd.DataFrame(overview)
    ov.to_csv(OUT / "Stage3A_cluster_marker_overview.csv", index=False)
    log("\nTOP CLUSTER MARKERS")
    log(ov.to_string(index=False))
log(f"Completed: {len(results)}/{len(SAMPLES)} samples")
log(
    "Interpretation: marker contrasts are within-specimen cluster-versus-rest descriptions, not validated cell identities."
)
log(
    "Pathology counts and candidate signatures require joint interpretation and spatial validation."
)
if errors:
    log("ERRORS:")
    [log(e) for e in errors]
else:
    log(
        "NEXT: review candidate markers and their spatial tissue context before pathway analysis."
    )
log(f"Outputs: {OUT}")
message = "\n".join(report)
(OUT / "Stage3A_Results.txt").write_text(message, encoding="utf-8")
print(message)
if errors:
    raise RuntimeError(f"Characterisation incomplete: {len(errors)} error(s); see Stage3B_Results.txt")

# ============================================================================
# Curated programme scoring and spatial summaries
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = PROJECT_ROOT
S1 = ROOT / "Results" / "Stage1D"
S2 = ROOT / "Results" / "Stage2A"
S2E = ROOT / "Results" / "Stage2E"
OUT = ROOT / "Results" / "Stage3B"
FIG = OUT / "Figures"
FIG.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

# Curated, interpretable gene sets; exploratory indicators, not cell-type calls.
PROGRAMMES = {
    "Epithelial": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT7",
        "TACSTD2",
        "MUC1",
        "CDH1",
    ],
    "T_cell_associated": [
        "CD3D",
        "CD3E",
        "CD3G",
        "TRAC",
        "CD2",
        "LCK",
        "IL7R",
        "CD247",
    ],
    "B_plasma_associated": [
        "MS4A1",
        "CD79A",
        "CD79B",
        "CD19",
        "MZB1",
        "JCHAIN",
        "IGKC",
        "IGHG1",
    ],
    "Stromal_ECM": [
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "DCN",
        "LUM",
        "FN1",
        "POSTN",
        "COL6A2",
    ],
    "Hypoxia_associated": [
        "NDUFA4L2",
        "NDRG1",
        "VEGFA",
        "SLC2A1",
        "BNIP3",
        "CA9",
        "EGLN3",
        "LDHA",
    ],
    "Proliferation_associated": [
        "MKI67",
        "TOP2A",
        "UBE2C",
        "CENPF",
        "PCNA",
        "TYMS",
        "MCM2",
        "STMN1",
    ],
}
MIN_PRESENT = 4
MIN_DETECTED_FRACTION = 0.01
MIN_SD = 0.03
report = []
summary = []
pathology_rows = []
cluster_rows = []
errors = []


def log(message=""):
    report.append(str(message))


def save_map(sample, label, obs, values):
    fig, ax = plt.subplots(figsize=(7.6, 6.9))
    x = pd.to_numeric(obs["array_col"]).to_numpy(dtype=float)
    y = pd.to_numeric(obs["array_row"]).to_numpy(dtype=float)
    finite = np.isfinite(values)
    low, high = np.percentile(values[finite], [2, 98])
    if low >= high:
        low, high = float(np.min(values)), float(np.max(values) + 1e-6)
    art = ax.scatter(
        x,
        y,
        c=values,
        s=14,
        marker="h",
        cmap="coolwarm",
        vmin=low,
        vmax=high,
        linewidths=0,
        rasterized=True,
    )
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(
        f"{sample} | {label.replace('_', ' ')}", fontweight="bold", fontsize=13
    )
    ax.set_xlabel("Visium array column", fontweight="bold")
    ax.set_ylabel("Visium array row", fontweight="bold")
    cb = fig.colorbar(art, ax=ax, shrink=0.78, pad=0.025)
    cb.set_label("Within-specimen mean gene z-score", fontweight="bold")
    fig.tight_layout()
    fig.savefig(
        FIG / f"{sample}_{label}_spatial_score.webp", dpi=250, bbox_inches="tight"
    )
    plt.close(fig)


def analyse(sample, genes):
    obs = pd.read_csv(S1 / f"{sample}_spot_metadata.csv", index_col="barcode")
    X = sparse.load_npz(S2 / f"{sample}_log1p_CP10K_spots_x_genes.npz").tocsr()
    assignment = pd.read_csv(S2E / f"{sample}_log_zscore_spot_clusters_UMAP.csv")
    if X.shape != (len(obs), len(genes)):
        raise ValueError(f"Expression shape mismatch: {X.shape}")
    if obs.index.has_duplicates or assignment["barcode"].duplicated().any():
        raise ValueError("Duplicate spot barcode")
    if not np.array_equal(
        obs.index.astype(str).to_numpy(), assignment["barcode"].astype(str).to_numpy()
    ):
        raise ValueError(
            "Spot order mismatched; not assigning scores to wrong coordinates"
        )
    if (X.data < 0).any() or not np.isfinite(X.data).all():
        raise ValueError("Nonfinite or negative log-normalised expression")
    for col in [
        "array_row",
        "array_col",
        "total_counts_shared_genes",
        "Classification",
    ]:
        if col not in obs:
            raise ValueError(f"Missing metadata column {col}")
    if not np.isfinite(obs[["array_row", "array_col"]].to_numpy(dtype=float)).all():
        raise ValueError("Nonfinite spatial coordinates")
    index = {g: i for i, g in enumerate(genes)}
    depth = np.log1p(obs["total_counts_shared_genes"].to_numpy(dtype=float))
    clusters = assignment["cluster"].to_numpy()
    pathology = obs["Classification"].fillna("[MISSING]").astype(str).to_numpy()
    scores = pd.DataFrame(
        {
            "barcode": obs.index,
            "cluster": clusters,
            "Classification": pathology,
            "array_row": obs["array_row"].to_numpy(),
            "array_col": obs["array_col"].to_numpy(),
        }
    )
    # Nearest spatial neighbours exclude self. Correlation is DESCRIPTIVE.
    xy = obs[["array_row", "array_col"]].to_numpy(dtype=float)
    knn = NearestNeighbors(n_neighbors=min(7, len(obs))).fit(xy)
    neighbors = knn.kneighbors(xy, return_distance=False)[:, 1:]
    for label, members in PROGRAMMES.items():
        present = [g for g in members if g in index]
        if not present:
            raise ValueError(f"{label}: no genes present; inspect gene annotation")
        a = X[:, [index[g] for g in present]].toarray().astype(np.float64)
        det = (a > 0).mean(axis=0)
        sd = a.std(axis=0, ddof=0)
        eligible = (det >= MIN_DETECTED_FRACTION) & (sd >= MIN_SD)
        used = [g for g, valid in zip(present, eligible) if valid]
        if len(used) < MIN_PRESENT:
            log(
                f"  WARNING {sample} {label}: {len(used)} usable genes; interpret cautiously"
            )
        if not eligible.any():
            raise ValueError(f"{label}: no variable, sufficiently detected genes")
        x = a[:, eligible]
        score = ((x - x.mean(axis=0)) / x.std(axis=0, ddof=0)).mean(axis=1)
        if not np.isfinite(score).all():
            raise ValueError("Nonfinite score")
        scores[label] = score.astype(np.float32)
        rdepth = float(spearmanr(score, depth).statistic)
        neighmean = score[neighbors].mean(axis=1)
        spatial_r = (
            float(np.corrcoef(score, neighmean)[0, 1])
            if np.std(score) > 0 and np.std(neighmean) > 0
            else np.nan
        )
        summary.append(
            dict(
                sample=sample,
                programme=label,
                genes_listed=len(members),
                genes_in_matrix=len(present),
                genes_used=len(used),
                used_genes=";".join(used),
                median_score=round(float(np.median(score)), 3),
                q25=round(float(np.percentile(score, 25)), 3),
                q75=round(float(np.percentile(score, 75)), 3),
                spearman_depth=round(rdepth, 3),
                spatial_neighbour_r=round(spatial_r, 3),
            )
        )
        for group in sorted(set(pathology)):
            m = pathology == group
            pathology_rows.append(
                dict(
                    sample=sample,
                    programme=label,
                    pathology=group,
                    n_spots=int(m.sum()),
                    median_score=round(float(np.median(score[m])), 3),
                    mean_score=round(float(np.mean(score[m])), 3),
                )
            )
        for group in sorted(set(clusters)):
            m = clusters == group
            cluster_rows.append(
                dict(
                    sample=sample,
                    programme=label,
                    cluster=group,
                    n_spots=int(m.sum()),
                    mean_score=round(float(np.mean(score[m])), 3),
                )
            )
        save_map(sample, label, obs, score)
    scores.to_csv(OUT / f"{sample}_spatial_programme_scores.csv", index=False)
    log(f"{sample}: PASS | spots={len(obs):,} | programme maps={len(PROGRAMMES)}")


log("PROJECT 2 | STAGE 3B — SPATIAL PROGRAMME MAPPING AND PATHOLOGY VALIDATION")
log(
    "Data: Stage2A log1p(CP10K), Stage1D coordinates/pathology, Stage2E original clusters."
)
log(
    "Score: mean within-specimen z-score across detected, variable curated marker genes."
)
log(
    "No inferred cell types, malignant status, spot-level p-values, or modifications of source data."
)
try:
    genes = (
        pd.read_csv(S1 / "shared_gene_symbols_ordered.csv")["gene_symbol"]
        .astype(str)
        .tolist()
    )
    if len(genes) != len(set(genes)):
        raise ValueError("Duplicate shared genes")
    for sample in SAMPLES:
        try:
            analyse(sample, genes)
        except Exception as e:
            errors.append(f"{sample}: {type(e).__name__}: {e}")
            log("ERROR: " + errors[-1])
            log(traceback.format_exc())
except Exception as e:
    errors.append(f"GLOBAL: {type(e).__name__}: {e}")
    log("ERROR: " + errors[-1])
    log(traceback.format_exc())

log("\n" + "=" * 92 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 92)
if summary:
    df = pd.DataFrame(summary)
    df.to_csv(OUT / "Stage3B_programme_summary.csv", index=False)
    log("\nPROGRAMME SCORES AND TECHNICAL ASSOCIATIONS")
    log(df.drop(columns="used_genes").to_string(index=False))
    log("\nGENES USED FOR EACH PROGRAMME")
    log(df[["sample", "programme", "used_genes"]].to_string(index=False))
if pathology_rows:
    p = pd.DataFrame(pathology_rows)
    p.to_csv(OUT / "Stage3B_pathology_programme_scores.csv", index=False)
    # Dominant and bottom pathology groups (min 10 spots) for concise interpretation.
    p10 = p[(p.n_spots >= 10) & (p.pathology != "[MISSING]")]
    log("\nPATHOLOGY: HIGHEST AND LOWEST MEDIAN SCORE (>=10 spots)")
    for (s, pr), g in p10.groupby(["sample", "programme"], sort=False):
        ordered = g.sort_values("median_score")
        lo = ordered.iloc[0]
        hi = ordered.iloc[-1]
        log(
            f"{s} | {pr}: LOW {lo.pathology} (n={lo.n_spots}, median={lo.median_score:+.3f}); "
            f"HIGH {hi.pathology} (n={hi.n_spots}, median={hi.median_score:+.3f})"
        )
if cluster_rows:
    pd.DataFrame(cluster_rows).to_csv(
        OUT / "Stage3B_cluster_programme_scores.csv", index=False
    )
log(f"Completed specimens: {len(set(r['sample'] for r in summary))}/{len(SAMPLES)}")
log(f"Spatial WebP maps saved: {len(summary)}")
log(
    "Caution: scores are relative WITHIN each tumour; do not compare absolute scores across tumours."
)
log(
    "Spatial neighbour correlations are descriptive, not spatial permutation-test p-values."
)
log(
    "Pathology descriptions are spot-level associations; tumours (n=6), not spots, are independent specimens."
)
if errors:
    log("ERRORS:")
    [log(e) for e in errors]
else:
    log(
        "NEXT: inspect spatial maps alongside histology, then investigate focused tumour–microenvironment hypotheses."
    )
log(f"Outputs: {OUT}")
message = "\n".join(report)
(OUT / "Stage3B_Results.txt").write_text(message, encoding="utf-8")
print(message)
