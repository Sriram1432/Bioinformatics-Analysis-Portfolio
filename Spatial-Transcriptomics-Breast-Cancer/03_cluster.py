"""Construct transcriptional neighbourhoods, UMAPs and Louvain clusters.

The primary analysis uses the established log-z-score PCA representation and
retains the alternative square-root-centred representation for sensitivity
comparison. Cluster labels remain specimen-specific descriptive labels.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ============================================================================
# Neighbour graphs, Louvain communities and UMAP coordinates
# ============================================================================

from pathlib import Path
import traceback
import warnings
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.sparse import csr_matrix

try:
    import umap
    import networkx as nx
except ImportError as err:
    raise SystemExit(
        f"Missing package: {err}. Install the required dependencies with "
        "'python -m pip install umap-learn networkx' and retry."
    )
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = PROJECT_ROOT
DATA = ROOT / "Results" / "Stage1D"
OLD = ROOT / "Results" / "Stage2B"
NEW = ROOT / "Results" / "Stage2D"
OUT = ROOT / "Results" / "Stage2E"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
SEED = 20261101
K = 15
PCS = 15
RESOLUTION = 0.7  # exploratory setting; no tuning to pathology labels
REPORT = []
SUMMARY = []
ERRORS = []
COMPARISONS = []


def log(s=""):
    REPORT.append(str(s))


def nmi_path(labels, pathology):
    p = pd.Series(pathology).fillna("[MISSING]").astype(str).to_numpy()
    mask = p != "[MISSING]"
    if len(np.unique(p[mask])) < 2:
        return np.nan
    return float(normalized_mutual_info_score(p[mask], labels[mask]))


def eta_depth(labels, depth):
    mu = depth.mean()
    ss = np.square(depth - mu).sum()
    return (
        float(
            sum(
                np.sum(labels == c) * np.square(depth[labels == c].mean() - mu)
                for c in np.unique(labels)
            )
            / ss
        )
        if ss > 0
        else np.nan
    )


def make_graph(z):
    nn = NearestNeighbors(n_neighbors=K + 1, metric="euclidean", n_jobs=-1)
    nn.fit(z)
    neighbor = nn.kneighbors(z, return_distance=False)
    rows = np.repeat(np.arange(len(z)), K)
    cols = neighbor[:, 1 : K + 1].ravel()
    # Unit weights: distances in the two PCA spaces have different scales.
    adjacency = csr_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(len(z), len(z))
    )
    adjacency = adjacency.maximum(adjacency.T)
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    graph = nx.from_scipy_sparse_array(adjacency)
    clusters = list(
        nx.algorithms.community.louvain_communities(
            graph, weight=None, resolution=RESOLUTION, seed=SEED
        )
    )
    labels = np.full(len(z), -1, dtype=int)
    for j, members in enumerate(sorted(clusters, key=lambda s: -len(s))):
        labels[list(members)] = j
    if np.any(labels < 0):
        raise ValueError("Unassigned clustering nodes")
    return labels


def spatial_coherence(obs, labels):
    # Spatial-neighbourhood test is DESCRIPTIVE only; tissue spot density is irregular.
    xy = obs[["array_row", "array_col"]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():
        raise ValueError("Nonfinite array coordinates")
    neighbors = NearestNeighbors(
        n_neighbors=min(7, len(labels)), metric="euclidean"
    ).fit(xy)
    ids = neighbors.kneighbors(xy, return_distance=False)[:, 1:]
    observed = float(np.mean(labels[:, None] == labels[ids]))
    frequencies = np.bincount(labels) / len(labels)
    random_baseline = float(np.sum(frequencies**2))
    return observed, random_baseline


def figure(sample, method, emb, labels, obs):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    axes[0].scatter(
        emb[:, 0], emb[:, 1], c=labels, s=3, cmap="tab20", alpha=0.8, rasterized=True
    )
    axes[0].set(
        title=f"{sample} | {method} | UMAP / Louvain clusters",
        xlabel="UMAP 1",
        ylabel="UMAP 2",
    )
    # Spatial cluster map (array coordinates, NOT UMAP coordinates).
    axes[1].scatter(
        obs["array_col"],
        obs["array_row"],
        c=labels,
        s=8,
        cmap="tab20",
        alpha=0.9,
        rasterized=True,
    )
    axes[1].invert_yaxis()
    axes[1].set_aspect("equal", adjustable="box")
    axes[1].set(
        title="Spatial cluster distribution (array coordinates)",
        xlabel="Array column",
        ylabel="Array row",
    )
    fig.savefig(OUT / f"{sample}_{method}_UMAP_Spatial_clusters.webp", dpi=200)
    plt.close(fig)


log("PROJECT 2 | STAGE 2E — UMAP + GRAPH CLUSTERING, TWO PCA REPRESENTATIONS")
log(
    f"Within each specimen: {PCS} PCs, {K} neighbours, Louvain resolution {RESOLUTION}, UMAP seed {SEED}."
)
log(
    "Comparisons are descriptive. No cell-type labels, inferential p-values, or cross-specimen spot pseudoreplication."
)
log(
    "UMAP = transcriptional neighbourhood view; spatial plots use physical array coordinates."
)

for sample in SAMPLES:
    try:
        obs = pd.read_csv(DATA / f"{sample}_spot_metadata.csv", index_col="barcode")
        baseline = pd.read_csv(
            OLD / f"{sample}_PCA30_spot_scores.csv", index_col="barcode"
        )
        alternative = pd.read_csv(
            NEW / f"{sample}_sqrtCP10K_PCA15_scores.csv", index_col="barcode"
        )
        if not (
            obs.index.equals(baseline.index) and obs.index.equals(alternative.index)
        ):
            raise ValueError("Barcode ordering mismatch: refusing to proceed")
        total = np.asarray(obs["total_counts_shared_genes"], float)
        if (total <= 0).any():
            raise ValueError("Nonpositive library sizes")
        depth = np.log1p(total)
        methods = {
            "log_zscore": baseline[[f"PC{i}" for i in range(1, PCS + 1)]].to_numpy(
                float
            ),
            "sqrt_centered": alternative[
                [f"PC{i}" for i in range(1, PCS + 1)]
            ].to_numpy(float),
        }
        assigned = {}
        for method, z in methods.items():
            if not np.isfinite(z).all():
                raise ValueError(f"{method}: nonfinite PCA scores")
            labels = make_graph(z)
            assigned[method] = labels
            embedding = umap.UMAP(
                n_neighbors=K,
                n_components=2,
                min_dist=0.3,
                metric="euclidean",
                random_state=SEED,
                n_jobs=1,
            ).fit_transform(z)
            if not np.isfinite(embedding).all():
                raise ValueError("Nonfinite UMAP output")
            counts = np.bincount(labels)
            p_nmi = nmi_path(labels, obs["Classification"].to_numpy())
            depth_eta = eta_depth(labels, depth)
            spatial_agree, spatial_null = spatial_coherence(obs, labels)
            out = pd.DataFrame(
                {
                    "barcode": obs.index,
                    "cluster": labels,
                    "UMAP1": embedding[:, 0],
                    "UMAP2": embedding[:, 1],
                }
            )
            out.to_csv(OUT / f"{sample}_{method}_spot_clusters_UMAP.csv", index=False)
            figure(sample, method, embedding, labels, obs)
            SUMMARY.append(
                dict(
                    sample=sample,
                    method=method,
                    spots=len(obs),
                    clusters=len(counts),
                    min_cluster_size=int(counts.min()),
                    max_cluster_size=int(counts.max()),
                    pathology_NMI=round(p_nmi, 3),
                    depth_eta2=round(depth_eta, 3),
                    spatial_neighbor_same=round(spatial_agree, 3),
                    spatial_random_baseline=round(spatial_null, 3),
                )
            )
        ari = adjusted_rand_score(assigned["log_zscore"], assigned["sqrt_centered"])
        COMPARISONS.append(
            dict(sample=sample, ARI_between_methods=round(float(ari), 3))
        )
        log(
            f"{sample}: PASS | spots={len(obs):,} | old clusters={len(np.unique(assigned['log_zscore']))} | "
            f"sqrt clusters={len(np.unique(assigned['sqrt_centered']))} | between-method ARI={ari:.3f}"
        )
    except Exception as exc:
        ERRORS.append(f"{sample}: {type(exc).__name__}: {exc}")
        log("ERROR: " + ERRORS[-1])
        log(traceback.format_exc())

log("\n" + "=" * 95 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 95)
if SUMMARY:
    table = pd.DataFrame(SUMMARY)
    table.to_csv(OUT / "Stage2E_clustering_diagnostics.csv", index=False)
    log(table.to_string(index=False))
if COMPARISONS:
    comp = pd.DataFrame(COMPARISONS)
    comp.to_csv(OUT / "Stage2E_method_agreement.csv", index=False)
    log("\nBETWEEN-METHOD CLUSTER AGREEMENT (ARI)")
    log(comp.to_string(index=False))
log(f"Completed: {len(COMPARISONS)}/{len(SAMPLES)} specimens")
log(
    "Caution: pathology NMI and depth eta² are DESCRIPTIVE; high values do not alone establish biology or confounding."
)
log(
    "Spatial-neighbour agreement should be compared with its cluster-frequency random baseline; neither is a spatial permutation test."
)
log(
    "Different representations may give different resolutions of one biological continuum. No winner selected automatically."
)
log("Existing files from Stages 1D–2D were not modified.")
if ERRORS:
    log("ERRORS:")
    for e in ERRORS:
        log(e)
else:
    log(
        "NEXT: review the 12 WebP figures and diagnostics, then choose a representation for marker programmes."
    )
log(f"Outputs: {OUT}")
result = "\n".join(REPORT)
(OUT / "Stage2E_Results.txt").write_text(result, encoding="utf-8")
print(result)
if ERRORS:
    raise RuntimeError(f"Clustering incomplete: {len(ERRORS)} error(s); see Stage2E_Results.txt")
