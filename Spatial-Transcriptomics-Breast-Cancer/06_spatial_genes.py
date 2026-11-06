"""Identify recurrent spatially structured genes within each tumour section.

Moran's I is calculated on a symmetric six-neighbour spatial graph for the 800
variance-ranked genes in each specimen. Empirical probabilities use fixed
permutations and Benjamini-Hochberg correction within specimen.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ============================================================================
# Moran's I calculation, permutation testing and recurrence summaries
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


ROOT = PROJECT_ROOT
S1 = ROOT / "Results" / "Stage1D"
S2 = ROOT / "Results" / "Stage2A"
S2B = ROOT / "Results" / "Stage2B"
OUT = ROOT / "Results" / "Stage5A"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
GENES_TESTED = 800  # Preselected by expression variance, NOT by spatial Moran I
PERMUTATIONS = 199  # Minimum empirical p-value = 1/200 = 0.005
NEIGHBOURS = 6
SEED = 20261101
BATCH = 160
report, summaries, failures = [], [], []


def log(s=""):
    report.append(str(s))


def bh(p):
    p = np.asarray(p, dtype=np.float64)
    o = np.argsort(p, kind="stable")
    q = np.empty(len(p), dtype=np.float64)
    q[o] = np.minimum(
        1, np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    )
    return q


def graph_from_coordinates(obs):
    xy = obs[["array_row", "array_col"]].to_numpy(dtype=np.float64)
    if not np.isfinite(xy).all() or np.unique(xy, axis=0).shape[0] != len(xy):
        raise ValueError("Missing or duplicate array coordinates")
    nn = NearestNeighbors(n_neighbors=NEIGHBOURS + 1).fit(xy)
    nearest = nn.kneighbors(xy, return_distance=False)[:, 1:]
    r = np.repeat(np.arange(len(xy)), NEIGHBOURS)
    c = nearest.ravel()
    W = sparse.csr_matrix(
        (np.ones(len(r), dtype=np.float32), (r, c)), shape=(len(xy), len(xy))
    )
    W = W.maximum(W.T).tocsr()
    W.setdiag(0)
    W.eliminate_zeros()
    W.data[:] = 1
    return W


def moran_for_matrix(W, A, factor):
    # Columns are already centred. Each column gets its own variance divisor.
    numerator = np.sum(A * (W @ A), axis=0, dtype=np.float64)
    denominator = np.sum(A * A, axis=0, dtype=np.float64)
    return factor * numerator / denominator


log("PROJECT 2 | STAGE 5A — SPATIALLY VARIABLE GENE DISCOVERY")
log(
    f"Within-tumour Moran I: top {GENES_TESTED} variance-ranked genes/sample; {NEIGHBOURS}-NN symmetric spatial graph."
)
log(
    f"{PERMUTATIONS} independent spot-label permutations per gene; BH correction within each specimen."
)
log(
    "Genes were prespecified from Stage2B expression variance, not selected by spatial I."
)
log(
    "Random spot-label permutations test spatial randomness, not depth-adjusted or autocorrelation-preserving nulls."
)

try:
    genes = (
        pd.read_csv(S1 / "shared_gene_symbols_ordered.csv")["gene_symbol"]
        .astype(str)
        .to_numpy()
    )
    if len(set(genes)) != len(genes):
        raise ValueError("Duplicated shared gene symbols")
    index = {g: i for i, g in enumerate(genes)}
except Exception as e:
    raise RuntimeError(f"Cannot read shared gene annotation: {e}") from e

for sample in SAMPLES:
    try:
        obs = pd.read_csv(S1 / f"{sample}_spot_metadata.csv", index_col="barcode")
        X = sparse.load_npz(S2 / f"{sample}_log1p_CP10K_spots_x_genes.npz").tocsr()
        if X.shape != (len(obs), len(genes)) or obs.index.has_duplicates:
            raise ValueError("Spot/gene dimensions invalid or duplicate barcodes")
        if not np.isfinite(X.data).all() or np.any(X.data < 0):
            raise ValueError("Invalid expression values")
        hvgs = pd.read_csv(S2B / f"{sample}_top2000_variable_genes.csv")
        if not {"gene_symbol", "variance_log_expression"}.issubset(hvgs.columns):
            raise ValueError("Required HVG columns missing")
        hvgs = hvgs.sort_values(
            ["variance_log_expression", "gene_symbol"], ascending=[False, True]
        )
        requested = hvgs["gene_symbol"].astype(str).head(GENES_TESTED).to_numpy()
        if len(requested) != GENES_TESTED or len(set(requested)) != GENES_TESTED:
            raise ValueError("Cannot uniquely select requested candidate genes")
        selected = np.array([index[g] for g in requested], dtype=int)
        W = graph_from_coordinates(obs)
        factor = len(obs) / float(W.sum())
        totals = pd.to_numeric(
            obs["total_counts_shared_genes"], errors="raise"
        ).to_numpy(dtype=np.float64)
        depth = np.log1p(totals)
        depth = (depth - depth.mean()) / depth.std()
        rng = np.random.default_rng(SEED + SAMPLES.index(sample))
        rows = []
        for start in range(0, GENES_TESTED, BATCH):
            stop = min(start + BATCH, GENES_TESTED)
            # A modest dense block avoids creating the entire dense transcriptome.
            A = X[:, selected[start:stop]].toarray().astype(np.float32)
            det = (A > 0).mean(axis=0)
            A -= A.mean(axis=0, keepdims=True)
            var = np.sum(A * A, axis=0, dtype=np.float64)
            if np.any(var <= 0):
                raise ValueError("Zero variance selected candidate gene")
            I = moran_for_matrix(W, A, factor)
            exceed = np.zeros(A.shape[1], dtype=np.int32)
            # Reuse each permutation across genes, giving separate conditional tests.
            for p in range(PERMUTATIONS):
                ix = rng.permutation(len(obs))
                Ip = moran_for_matrix(W, A[ix, :], factor)
                exceed += Ip >= I - 1e-12
            pvals = (exceed + 1) / (PERMUTATIONS + 1)
            # Pearson correlation with log depth: technical sensitivity flag.
            depth_r = (depth.astype(np.float64) @ A.astype(np.float64)) / (
                np.sqrt(np.sum(A * A, axis=0, dtype=np.float64)) * np.linalg.norm(depth)
            )
            for j, g in enumerate(requested[start:stop]):
                rows.append(
                    {
                        "sample": sample,
                        "gene_symbol": g,
                        "moran_I": float(I[j]),
                        "empirical_p_positive": float(pvals[j]),
                        "fraction_spots_detected": float(det[j]),
                        "pearson_r_log_depth": float(depth_r[j]),
                    }
                )
        result = pd.DataFrame(rows)
        result["BH_FDR"] = bh(result["empirical_p_positive"])
        result["FDR_lt_0_05"] = (result.BH_FDR < 0.05) & (result.moran_I > 0)
        result = result.sort_values(["moran_I", "gene_symbol"], ascending=[False, True])
        result.to_csv(OUT / f"{sample}_MoranI_gene_results.csv", index=False)
        sig = result[result.FDR_lt_0_05]
        record = {
            "sample": sample,
            "spots": len(obs),
            "genes_tested": len(result),
            "spatial_graph_edges": W.nnz // 2,
            "FDR_significant_positive": len(sig),
            "median_moran_I": float(result.moran_I.median()),
            "max_moran_I": float(result.moran_I.max()),
            "significant_abs_depth_r_gt_0_5": int(
                (sig.pearson_r_log_depth.abs() > 0.5).sum()
            ),
        }
        summaries.append(record)
        log(
            f"{sample}: PASS | {len(obs):,} spots | {len(result)} tested | "
            f"{len(sig)} positive FDR<0.05 | median I={record['median_moran_I']:.3f}"
        )
        top = result.head(10)
        log(
            "  Top spatial genes: "
            + ", ".join(
                f"{r.gene_symbol} (I={r.moran_I:.3f}, q={r.BH_FDR:.3g})"
                for r in top.itertuples()
            )
        )
        del X, W
    except Exception as exc:
        failures.append(f"{sample}: {type(exc).__name__}: {exc}")
        log("ERROR: " + failures[-1])
        log(traceback.format_exc())

log("\n" + "=" * 86 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 86)
if summaries:
    tab = pd.DataFrame(summaries)
    tab.to_csv(OUT / "Stage5A_summary.csv", index=False)
    log(tab.to_string(index=False))
    if len(summaries) == len(SAMPLES):
        joined = pd.concat(
            [pd.read_csv(OUT / f"{s}_MoranI_gene_results.csv") for s in SAMPLES],
            ignore_index=True,
        )
        recurrent = (
            joined[joined.FDR_lt_0_05]
            .groupby("gene_symbol")
            .agg(
                specimens=("sample", "nunique"),
                median_I=("moran_I", "median"),
                max_abs_depth_r=(
                    "pearson_r_log_depth",
                    lambda v: float(np.max(np.abs(v))),
                ),
            )
            .reset_index()
            .sort_values(["specimens", "median_I"], ascending=False)
        )
        recurrent.to_csv(OUT / "Stage5A_recurrent_spatial_genes.csv", index=False)
        log(
            "\nRECURRENT POSITIVE SPATIAL GENES (conditional on variance-based gene selection)"
        )
        log(
            recurrent.head(25).to_string(index=False)
            if len(recurrent)
            else "None at FDR < 0.05"
        )
        log(
            f"Genes significant in >=3 tumours: {int((recurrent.specimens >= 3).sum())}"
        )
log(f"Completed specimens: {len(summaries)}/6")
log(
    "LIMITATIONS: tests only 800 expression-variable genes/tumour (NOT genome-wide); Monte Carlo p min=0.005."
)
log(
    "BH is across tested genes per tumour; p-values refer to spot-label randomisation, NOT independent patients."
)
log(
    "Spatial structure may reflect pathology, sequencing-depth gradients, histological composition or technical features."
)
log(
    "Cross-tumour recurrence is descriptive; candidate lists differ by within-tumour variance selection."
)
log(f"Outputs: {OUT}")
if failures:
    log("ERRORS:\n" + "\n".join(failures))
(OUT / "Stage5A_Results.txt").write_text("\n".join(report), encoding="utf-8")
print("\n".join(report))
if failures:
    raise RuntimeError(
        f"Stage 5A incomplete: {len(failures)} samples failed; see report."
    )
