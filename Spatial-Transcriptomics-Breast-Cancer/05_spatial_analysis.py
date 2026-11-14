"""Run pathology-anchored spatial analyses and focused sensitivity checks.

The analyses cover hypoxia around necrosis, immune-associated distance patterns,
DCIS programme contrasts and robustness to neighbourhood or sequencing-depth
choices. All conclusions remain descriptive and specimen-aware.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ============================================================================
# Hypoxia-associated expression around necrotic tissue
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = PROJECT_ROOT
SAMPLE = "1142243F"
META = ROOT / "Results" / "Stage1D" / f"{SAMPLE}_spot_metadata.csv"
SCORES = ROOT / "Results" / "Stage3B" / f"{SAMPLE}_spatial_programme_scores.csv"
OUT = ROOT / "Results" / "Stage4A"
FIG = OUT / "Figures"
FIG.mkdir(parents=True, exist_ok=True)
SEED = 20261101
RNG = np.random.default_rng(SEED)
N_PERM = 999
BLOCK_WIDTH = (
    8  # Visium array-grid units. Spatial sensitivity, not true spatial-null test.
)
RANGE_BINS = [0, 2, 5, 10, np.inf]
RANGE_LABELS = ["<2", "2–<5", "5–<10", ">=10"]
REPORT = []


def log(s=""):
    REPORT.append(str(s))


def save_scatter(df, variable, name, title, cmap="viridis"):
    fig, ax = plt.subplots(figsize=(8, 6.6))
    c = df[variable].to_numpy(dtype=float)
    if np.isfinite(c).any():
        lo, hi = np.percentile(c[np.isfinite(c)], [2, 98])
        if lo >= hi:
            lo, hi = float(np.nanmin(c)), float(np.nanmax(c) + 1e-6)
        collection = ax.scatter(
            df.array_col,
            df.array_row,
            c=c,
            s=13,
            cmap=cmap,
            vmin=lo,
            vmax=hi,
            linewidths=0,
            rasterized=True,
        )
        cb = fig.colorbar(collection, ax=ax, shrink=0.8)
        cb.set_label(variable.replace("_", " "), fontweight="bold")
    ax.scatter(
        df.loc[df.is_necrosis, "array_col"],
        df.loc[df.is_necrosis, "array_row"],
        facecolors="none",
        edgecolors="black",
        s=30,
        linewidths=0.65,
        label="Necrotic annotated spot",
    )
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Visium array column", fontweight="bold")
    ax.set_ylabel("Visium array row", fontweight="bold")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=260, bbox_inches="tight")
    plt.close(fig)


def block_permutation_test(df, xcol, ycol, nperm=N_PERM):
    """Spatially grouped score re-assignment as a sensitivity check only.

    Blocks are aggregated then their mean scores permuted across blocks;
    intentionally does NOT claim to preserve spatial covariance perfectly.
    """
    work = df[[xcol, ycol, "array_col", "array_row"]].dropna().copy()
    col = (work.array_col.astype(float) - work.array_col.min()) // BLOCK_WIDTH
    row = (work.array_row.astype(float) - work.array_row.min()) // BLOCK_WIDTH
    work["block"] = row.astype(int).astype(str) + "_" + col.astype(int).astype(str)
    blocks = work.groupby("block", sort=False).agg(
        distance=(xcol, "median"), score=(ycol, "median"), n=(ycol, "size")
    )
    if len(blocks) < 8:
        return np.nan, np.nan, len(blocks)
    observed = float(spearmanr(blocks.distance, blocks.score).statistic)
    if not np.isfinite(observed):
        return np.nan, np.nan, len(blocks)
    null = np.empty(nperm, dtype=float)
    for i in range(nperm):
        null[i] = spearmanr(
            blocks.distance, RNG.permutation(blocks.score.to_numpy())
        ).statistic
    p = (1 + np.count_nonzero(np.abs(null) >= abs(observed))) / (nperm + 1)
    return observed, float(p), len(blocks)


log("PROJECT 2 | STAGE 4A — HYPOXIA AND NECROTIC TISSUE SPATIAL GRADIENT")
log(
    "Sample: 1142243F | exploratory single-specimen study | no independent healthy controls"
)
log(
    "Distances: Visium array coordinates, not micrometres. Necrosis masks: pathology labels."
)
log("No changes to raw matrices, annotations or previously validated outputs.")
try:
    obs = pd.read_csv(META, index_col="barcode")
    score = pd.read_csv(SCORES)
    required = [
        "barcode",
        "Hypoxia_associated",
        "Epithelial",
        "Proliferation_associated",
    ]
    for c in required:
        if c not in score:
            raise ValueError(f"Stage3B scores missing column: {c}")
    for c in ["array_row", "array_col", "Classification", "total_counts_shared_genes"]:
        if c not in obs:
            raise ValueError(f"Stage1D metadata missing column: {c}")
    if obs.index.has_duplicates or score.barcode.duplicated().any():
        raise ValueError("Duplicate barcodes")
    if len(obs) != len(score) or set(obs.index.astype(str)) != set(
        score.barcode.astype(str)
    ):
        raise ValueError("Spot barcodes do not match across Stage1D and Stage3B")
    # Explicit key-based merge ensures no positional misalignment.
    df = obs.reset_index().merge(
        score[required], on="barcode", validate="one_to_one", how="left"
    )
    for col in [
        "array_row",
        "array_col",
        "Hypoxia_associated",
        "Epithelial",
        "Proliferation_associated",
        "total_counts_shared_genes",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if not np.isfinite(
        df[
            [
                "array_row",
                "array_col",
                "Hypoxia_associated",
                "Epithelial",
                "Proliferation_associated",
                "total_counts_shared_genes",
            ]
        ].to_numpy()
    ).all():
        raise ValueError("Missing/nonfinite coordinates, scores, or counts")
    if not (df.total_counts_shared_genes > 0).all():
        raise ValueError("Nonpositive spot depth")
    df["is_necrosis"] = (
        df.Classification.fillna("").str.strip().str.casefold().eq("necrosis")
    )
    xy = df[["array_row", "array_col"]].to_numpy(dtype=float)
    mask = df.is_necrosis.to_numpy()
    if mask.sum() < 10 or (~mask).sum() < 10:
        raise ValueError(
            "Insufficient necrosis or non-necrosis spots for distance analysis"
        )
    # For any non-necrotic spot: distance to nearest necrotic spot centre.
    # For necrotic spots: distance to closest non-necrotic spot centre (negative).
    nearest_necrosis = cKDTree(xy[mask])
    nearest_other = cKDTree(xy[~mask])
    d = np.empty(len(df), dtype=float)
    d[~mask] = nearest_necrosis.query(xy[~mask], k=1)[0]
    d[mask] = -nearest_other.query(xy[mask], k=1)[0]
    df["signed_necrosis_distance_array"] = d
    non = df.loc[~df.is_necrosis].copy()
    non["distance_bin"] = pd.cut(
        non.signed_necrosis_distance_array,
        bins=RANGE_BINS,
        labels=RANGE_LABELS,
        right=False,
        include_lowest=True,
    )
    df["log_depth"] = np.log1p(df.total_counts_shared_genes)
    non["log_depth"] = np.log1p(non.total_counts_shared_genes)
    # Note: these are selected summary endpoints, not 15,399 independent patients.
    cats = pd.concat(
        [df.loc[mask].assign(region="Necrosis"), non.assign(region="Non-necrosis")],
        ignore_index=True,
    )
    overall = (
        cats.groupby("region", sort=False)
        .agg(
            spots=("Hypoxia_associated", "size"),
            hypoxia_median=("Hypoxia_associated", "median"),
            epithelial_median=("Epithelial", "median"),
            proliferation_median=("Proliferation_associated", "median"),
        )
        .reset_index()
    )
    bands = (
        non.groupby("distance_bin", observed=False)
        .agg(
            spots=("Hypoxia_associated", "size"),
            hypoxia_median=("Hypoxia_associated", "median"),
            hypoxia_mean=("Hypoxia_associated", "mean"),
            depth_median=("total_counts_shared_genes", "median"),
            epithelial_median=("Epithelial", "median"),
            proliferation_median=("Proliferation_associated", "median"),
        )
        .reset_index()
    )
    if (bands.spots == 0).any():
        log("WARNING: one or more distance bins empty; interpret only populated bins.")
    rho = float(
        spearmanr(non.signed_necrosis_distance_array, non.Hypoxia_associated).statistic
    )
    rho_depth = float(
        spearmanr(non.signed_necrosis_distance_array, non.log_depth).statistic
    )
    rho_hyp_depth = float(spearmanr(non.Hypoxia_associated, non.log_depth).statistic)
    block_r, block_p, nblocks = block_permutation_test(
        non, "signed_necrosis_distance_array", "Hypoxia_associated"
    )
    # Depth-stratified DESCRIPTIVE correlations, to check gradient direction.
    non["depth_quartile"] = (
        pd.qcut(non.log_depth.rank(method="first"), 4, labels=False) + 1
    )
    depth_rows = []
    for q, g in non.groupby("depth_quartile"):
        depth_rows.append(
            dict(
                depth_quartile=int(q),
                spots=len(g),
                distance_hypoxia_spearman=round(
                    float(
                        spearmanr(
                            g.signed_necrosis_distance_array, g.Hypoxia_associated
                        ).statistic
                    ),
                    4,
                ),
            )
        )
    depth_table = pd.DataFrame(depth_rows)
    # Save full per-spot data and compact outputs.
    df.to_csv(OUT / f"{SAMPLE}_signed_distance_and_scores.csv", index=False)
    overall.to_csv(OUT / "Stage4A_necrosis_vs_other.csv", index=False)
    bands.to_csv(OUT / "Stage4A_distance_bins.csv", index=False)
    depth_table.to_csv(OUT / "Stage4A_depth_sensitivity.csv", index=False)
    save_scatter(
        df,
        "Hypoxia_associated",
        f"{SAMPLE}_hypoxia_necrosis.webp",
        f"{SAMPLE} | Hypoxia signature and necrosis",
        "coolwarm",
    )
    save_scatter(
        df,
        "signed_necrosis_distance_array",
        f"{SAMPLE}_necrosis_distance.webp",
        f"{SAMPLE} | Signed distance to necrosis",
        "viridis",
    )
    fig, ax = plt.subplots(figsize=(7.8, 5.5))
    valid = bands[bands.spots > 0]
    ax.plot(np.arange(len(valid)), valid.hypoxia_median, marker="o", linewidth=2)
    ax.set_xticks(np.arange(len(valid)), valid.distance_bin.astype(str))
    ax.set_xlabel("Distance to nearest necrotic spot (array units)", fontweight="bold")
    ax.set_ylabel("Median within-tumour hypoxia signature score", fontweight="bold")
    ax.set_title(f"{SAMPLE} | Non-necrotic tissue distance gradient", fontweight="bold")
    for i, (_, r) in enumerate(valid.iterrows()):
        ax.annotate(
            f"n={int(r.spots)}",
            (i, r.hypoxia_median),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=9,
        )
    fig.tight_layout()
    fig.savefig(FIG / f"{SAMPLE}_distance_gradient.webp", dpi=260, bbox_inches="tight")
    plt.close(fig)
    log(
        f"Matched {len(df):,} spots; necrosis={int(mask.sum()):,}; other={int((~mask).sum()):,}"
    )
    log("\nNECROSIS VERSUS NON-NECROSIS — DESCRIPTIVE MEDIANS")
    log(overall.round(4).to_string(index=False))
    log("\nNON-NECROTIC SPOTS — DISTANCE TO NEAREST NECROTIC SPOT")
    log(bands.round(4).to_string(index=False))
    log("\nPRIMARY SPOT-LEVEL ASSOCIATIONS (descriptive, spatially dependent)")
    log(f"Spearman hypoxia vs distance, NON-necrotic: {rho:+.4f}")
    log(f"Spearman depth vs distance, NON-necrotic: {rho_depth:+.4f}")
    log(f"Spearman hypoxia vs log depth, NON-necrotic: {rho_hyp_depth:+.4f}")
    log("\nDEPTH QUARTILE SENSITIVITY")
    log(depth_table.to_string(index=False))
    log(
        "\nEXPLORATORY SPATIAL-BLOCK SENSITIVITY (block medians; block locations shuffled)"
    )
    log(
        f"Block width={BLOCK_WIDTH} array units; blocks={nblocks}; permutations={N_PERM}"
    )
    log(f"Block-median Spearman rho={block_r:+.4f}; shuffled-block p={block_p:.4f}")
    log("WARNING: block-shuffle p is NOT a spatial autocorrelation-preserving null;")
    log(
        "it is a diagnostic, not definitive spatial significance or patient-level inference."
    )
    log(
        "Array grid distance approximates proximity of spot centres; it does not measure true necrotic boundary."
    )
    log(
        "Histology review required: some hypoxia scores may arise from necrosis-related tissue quality."
    )
    log(
        "Biological replication is one tumour here; do not generalise to patients from spot-level counts."
    )
    log(f"Figures saved: 3 WebP | Outputs: {OUT}")
    log(
        "NEXT: use actual result direction and spatial maps to choose immune–tumour investigation."
    )
except Exception as e:
    log(f"ERROR: {type(e).__name__}: {e}")
    log(traceback.format_exc())
finally:
    result = "\n".join(REPORT)
    (OUT / "Stage4A_Results.txt").write_text(result, encoding="utf-8")
    print(
        "\n" + "=" * 80 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 80
    )
    print(result)

# ============================================================================
# Immune-associated expression around lymphocyte annotations
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = PROJECT_ROOT
SAMPLES = ["1160920F", "CID44971", "CID4535"]
INPUT_META = ROOT / "Results" / "Stage1D"
INPUT_SCORES = ROOT / "Results" / "Stage3B"
OUT = ROOT / "Results" / "Stage4B"
FIG = OUT / "Figures"
FIG.mkdir(parents=True, exist_ok=True)
SIGNATURES = ["T_cell_associated", "B_plasma_associated", "Epithelial", "Stromal_ECM"]
BINS = [0, 2, 5, 10, np.inf]
LABELS = ["<2", "2–<5", "5–<10", ">=10"]
REPORT = []
TABLES = []
GROUPS = []
ERRORS = []


def log(s=""):
    REPORT.append(str(s))


def rank_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    good = np.isfinite(a) & np.isfinite(b)
    if good.sum() < 5 or np.ptp(a[good]) == 0 or np.ptp(b[good]) == 0:
        return np.nan
    return float(spearmanr(a[good], b[good]).statistic)


def plot_signature(df, sample, variable, label):
    x, y = df["array_col"].to_numpy(), df["array_row"].to_numpy()
    v = df[variable].to_numpy(dtype=float)
    lo, hi = np.percentile(v, [2, 98])
    if hi <= lo:
        hi = lo + 1e-6
    fig, ax = plt.subplots(figsize=(8, 6.7))
    im = ax.scatter(
        x,
        y,
        c=v,
        s=15,
        marker="h",
        cmap="coolwarm",
        vmin=lo,
        vmax=hi,
        edgecolors="none",
        rasterized=True,
    )
    # Black rings show pathology-annotated, lymphocyte-only anchor spots.
    m = df["lymphocyte_anchor"].to_numpy()
    ax.scatter(
        x[m],
        y[m],
        s=37,
        facecolors="none",
        edgecolors="black",
        linewidths=0.7,
        label="Lymphocytes (pathology)",
    )
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Visium array column", fontweight="bold")
    ax.set_ylabel("Visium array row", fontweight="bold")
    ax.set_title(f"{sample} | {label}", fontweight="bold")
    ax.legend(loc="best", fontsize=8, frameon=False)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("Within-specimen signature score", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / f"{sample}_{variable}_spatial.webp", dpi=250, bbox_inches="tight")
    plt.close(fig)


def plot_gradient(band_df, sample):
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    for signature in ["T_cell_associated", "B_plasma_associated", "Epithelial"]:
        s = band_df[band_df["programme"] == signature].set_index("distance_bin")
        v = [s.loc[b, "median_score"] if b in s.index else np.nan for b in LABELS]
        ax.plot(
            np.arange(len(LABELS)),
            v,
            marker="o",
            linewidth=2,
            label=signature.replace("_", " "),
        )
    ax.set_xticks(range(len(LABELS)), LABELS)
    ax.set_xlabel(
        "Distance to nearest lymphocyte-only spot (array units)", fontweight="bold"
    )
    ax.set_ylabel("Median within-tumour score", fontweight="bold")
    ax.set_title(f"{sample} | Non-anchor distance profiles", fontweight="bold")
    ax.legend(fontsize=9)
    ax.axhline(0, linewidth=0.8, color="gray")
    fig.tight_layout()
    fig.savefig(
        FIG / f"{sample}_immune_distance_profiles.webp", dpi=250, bbox_inches="tight"
    )
    plt.close(fig)


def analyse(sample):
    obs = pd.read_csv(INPUT_META / f"{sample}_spot_metadata.csv", index_col="barcode")
    score = pd.read_csv(INPUT_SCORES / f"{sample}_spatial_programme_scores.csv")
    if obs.index.has_duplicates or score["barcode"].duplicated().any():
        raise ValueError("Duplicate barcodes")
    if not set(SIGNATURES).issubset(score.columns):
        raise ValueError(
            f"Missing Stage3B scores: {set(SIGNATURES) - set(score.columns)}"
        )
    if "Classification" not in obs or not {
        "array_row",
        "array_col",
        "total_counts_shared_genes",
    }.issubset(obs.columns):
        raise ValueError("Stage1D metadata missing expected columns")
    if len(obs) != len(score) or set(obs.index.astype(str)) != set(
        score["barcode"].astype(str)
    ):
        raise ValueError("Stage1D and Stage3B barcode mismatch")
    df = obs.reset_index().merge(
        score[["barcode"] + SIGNATURES], how="left", on="barcode", validate="one_to_one"
    )
    required = ["array_row", "array_col", "total_counts_shared_genes"] + SIGNATURES
    for c in required:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if not np.isfinite(df[required].to_numpy(dtype=float)).all():
        raise ValueError(
            "Missing/nonfinite coordinates, library sizes or signature scores"
        )
    if not (df["total_counts_shared_genes"] > 0).all():
        raise ValueError("Nonpositive library size")
    # Only exact 'Lymphocytes' labels are anchor sites. Mixed categories are not
    # pure immune or pure tumour; we retain them as separate descriptive groups.
    labels = df["Classification"].fillna("[MISSING]").astype(str).str.strip()
    df["lymphocyte_anchor"] = labels.str.casefold().eq("lymphocytes")
    n_anchor = int(df["lymphocyte_anchor"].sum())
    if n_anchor < 10:
        raise ValueError(
            f"Only {n_anchor} lymphocyte-only labelled spots; cannot anchor distances"
        )
    xy = df[["array_row", "array_col"]].to_numpy(dtype=float)
    tree = cKDTree(xy[df["lymphocyte_anchor"].to_numpy()])
    dist = tree.query(xy, k=1)[0]
    df["distance_to_lymphocytes_array_units"] = dist
    non = df.loc[~df["lymphocyte_anchor"]].copy()
    non["distance_bin"] = pd.cut(
        non["distance_to_lymphocytes_array_units"],
        bins=BINS,
        labels=LABELS,
        include_lowest=True,
        right=False,
    )
    df["pathology"] = labels.to_numpy()
    grouped = []
    for category, subset in df.groupby("pathology", sort=True):
        row = {"sample": sample, "pathology": category, "spots": len(subset)}
        for signature in SIGNATURES:
            row[signature + "_median"] = round(float(subset[signature].median()), 4)
        grouped.append(row)
    GROUPS.extend(grouped)
    for signature in SIGNATURES:
        s = non[signature].to_numpy(dtype=float)
        dist_non = non["distance_to_lymphocytes_array_units"].to_numpy(dtype=float)
        depth = np.log1p(non["total_counts_shared_genes"].to_numpy(dtype=float))
        corr_dist = rank_corr(dist_non, s)
        corr_depth = rank_corr(depth, s)
        groups = non.groupby("distance_bin", observed=False)
        for band in LABELS:
            subset = (
                groups.get_group(band)
                if band in groups.groups and len(groups.groups[band])
                else non.iloc[:0]
            )
            TABLES.append(
                dict(
                    sample=sample,
                    programme=signature,
                    distance_bin=band,
                    spots=int(len(subset)),
                    median_score=round(float(subset[signature].median()), 4)
                    if len(subset)
                    else np.nan,
                    mean_score=round(float(subset[signature].mean()), 4)
                    if len(subset)
                    else np.nan,
                    median_log_depth=round(
                        float(np.log1p(subset["total_counts_shared_genes"]).median()), 4
                    )
                    if len(subset)
                    else np.nan,
                )
            )
        anchor_median = float(df.loc[df.lymphocyte_anchor, signature].median())
        far = non[non["distance_to_lymphocytes_array_units"] >= 10]
        far_median = float(far[signature].median()) if len(far) else np.nan
        log(
            f"  {signature}: anchor median={anchor_median:+.3f}; distant >=10 median={far_median:+.3f};"
            f" rho(distance,score)={corr_dist:+.3f}; rho(log-depth,score)={corr_depth:+.3f}"
        )
    for signature in ["T_cell_associated", "B_plasma_associated"]:
        plot_signature(df, sample, signature, signature.replace("_", " "))
    bands = pd.DataFrame([x for x in TABLES if x["sample"] == sample])
    plot_gradient(bands, sample)
    df.to_csv(OUT / f"{sample}_immune_distances_and_scores.csv", index=False)
    log(
        f"{sample}: PASS | spots={len(df):,} | lymphocyte-only anchors={n_anchor} | non-anchor={len(non):,}"
    )
    log(
        "  Pathology categories: "
        + "; ".join(f"{k}={v}" for k, v in labels.value_counts().items())
    )


log("PROJECT 2 | STAGE 4B — TUMOUR MICROENVIRONMENT AND IMMUNE-SPATIAL RELATIONSHIPS")
log(
    'Exact pathology label "Lymphocytes" is used as the spatial anchor; all other labels are retained descriptively.'
)
log(
    "Distances are Visium array grid units (not micrometres); spot labels do not establish individual cell identity."
)
log(
    "No inferential p-values: one patient per specimen, spatial dependence among spots."
)
log(
    "Signatures are relative within each specimen and should not be compared as absolute values between patients."
)
for sample in SAMPLES:
    try:
        analyse(sample)
    except Exception as exc:
        ERRORS.append(f"{sample}: {type(exc).__name__}: {exc}")
        log("ERROR: " + ERRORS[-1])
        log(traceback.format_exc())
log("\n" + "=" * 85 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 85)
if TABLES:
    bands = pd.DataFrame(TABLES)
    bands.to_csv(OUT / "Stage4B_immune_distance_bands.csv", index=False)
    log("\nDISTANCE BANDS — NON-ANCHOR SPOTS")
    log(
        bands[
            ["sample", "programme", "distance_bin", "spots", "median_score"]
        ].to_string(index=False)
    )
if GROUPS:
    pathology = pd.DataFrame(GROUPS)
    pathology.to_csv(OUT / "Stage4B_pathology_signature_medians.csv", index=False)
    log("\nPATHOLOGY SIGNATURE MEDIANS (at least 10 spots)")
    log(pathology[pathology.spots >= 10].to_string(index=False))
log(f"\nCompleted specimens: {len(SAMPLES) - len(ERRORS)}/{len(SAMPLES)}")
log(
    f"Saved spatial WebP figures: {3 * (len(SAMPLES) - len(ERRORS))} if all plot operations succeeded"
)
log(
    "Caution: correlations are descriptive, not patient-level evidence or spatial-autocorrelation-adjusted tests."
)
log(
    "High signature scores in mixed spots do not establish immune infiltration into malignant cells."
)
if ERRORS:
    log("\nERRORS:\n" + "\n".join(ERRORS))
log(f"Outputs: {OUT}")
text = "\n".join(REPORT)
(OUT / "Stage4B_Consolidated_Results.txt").write_text(text, encoding="utf-8")
print(text)

# ============================================================================
# DCIS and mixed invasive-associated programme contrasts
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = PROJECT_ROOT
SAMPLE = "CID44971"
OUT = ROOT / "Results" / "Stage4C"
FIG = OUT / "Figures"
FIG.mkdir(parents=True, exist_ok=True)
PROGRAMMES = [
    "Epithelial",
    "Proliferation_associated",
    "Hypoxia_associated",
    "T_cell_associated",
    "B_plasma_associated",
    "Stromal_ECM",
]
GROUPS = {"DCIS": "DCIS", "Invasive cancer + lymphocytes": "Invasive-associated mixed"}
REPORT = []


def log(s=""):
    REPORT.append(str(s))


def med(x):
    a = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if len(a) else np.nan


def corr(a, b):
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 5 or np.ptp(x[valid]) == 0 or np.ptp(y[valid]) == 0:
        return np.nan
    return float(spearmanr(x[valid], y[valid]).statistic)


try:
    meta_path = ROOT / "Results" / "Stage1D" / f"{SAMPLE}_spot_metadata.csv"
    score_path = ROOT / "Results" / "Stage3B" / f"{SAMPLE}_spatial_programme_scores.csv"
    meta = pd.read_csv(meta_path, index_col="barcode")
    scores = pd.read_csv(score_path)
    if meta.index.has_duplicates or scores["barcode"].duplicated().any():
        raise ValueError("Duplicate spot barcodes detected")
    needed_meta = [
        "Classification",
        "array_col",
        "array_row",
        "total_counts_shared_genes",
    ]
    for col in needed_meta:
        if col not in meta:
            raise ValueError(f"Missing metadata field: {col}")
    for col in PROGRAMMES:
        if col not in scores:
            raise ValueError(f"Missing Stage3B score: {col}")
    if len(meta) != len(scores) or set(meta.index.astype(str)) != set(
        scores.barcode.astype(str)
    ):
        raise ValueError("Stage1D and Stage3B spot barcodes do not match exactly")
    df = meta.reset_index().merge(
        scores[["barcode"] + PROGRAMMES], on="barcode", validate="one_to_one"
    )
    for col in ["array_col", "array_row", "total_counts_shared_genes"] + PROGRAMMES:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if not np.isfinite(
        df[
            ["array_col", "array_row", "total_counts_shared_genes"] + PROGRAMMES
        ].to_numpy(dtype=float)
    ).all():
        raise ValueError("Missing or nonfinite coordinates, depth or scores")
    if (df.total_counts_shared_genes <= 0).any():
        raise ValueError("Nonpositive library size")
    labels = df.Classification.fillna("[MISSING]").astype(str).str.strip()
    df["pathology"] = labels
    focus = df[df.pathology.isin(GROUPS)].copy()
    focus["group"] = focus.pathology.map(GROUPS)
    if set(focus.pathology.unique()) != set(GROUPS):
        raise ValueError(
            "Expected DCIS and invasive-associated pathology groups not both present"
        )
    counts = focus.groupby("group").size().to_dict()
    log("PROJECT 2 | STAGE 4C — DCIS VS INVASIVE-ASSOCIATED TISSUE PROGRAMMES")
    log(
        "One patient (CID44971). Comparisons are DESCRIPTIVE; spots are spatially dependent."
    )
    log(
        'DCIS versus "Invasive cancer + lymphocytes": groups differ in cell composition.'
    )
    log(
        "Differences CANNOT establish invasion-specific effects or malignant-cell states."
    )
    log(
        "Scores are within-specimen z-score averages; differences are in score units, not log fold changes."
    )
    log(
        f"Input verified: {len(df):,} spots | DCIS={counts.get('DCIS', 0)} | invasive-associated mixed={counts.get('Invasive-associated mixed', 0)}"
    )
    rows = []
    for p in PROGRAMMES:
        a = focus.loc[focus.group == "DCIS", p]
        b = focus.loc[focus.group == "Invasive-associated mixed", p]
        rows.append(
            {
                "programme": p,
                "DCIS_n": len(a),
                "mixed_invasive_n": len(b),
                "DCIS_median": round(med(a), 4),
                "mixed_invasive_median": round(med(b), 4),
                "difference_DCIS_minus_mixed": round(med(a) - med(b), 4),
                "DCIS_Q1": round(float(a.quantile(0.25)), 4),
                "DCIS_Q3": round(float(a.quantile(0.75)), 4),
                "mixed_Q1": round(float(b.quantile(0.25)), 4),
                "mixed_Q3": round(float(b.quantile(0.75)), 4),
                "rho_with_log_depth_all_spots": round(
                    corr(np.log1p(df.total_counts_shared_genes), df[p]), 4
                ),
            }
        )
    results = pd.DataFrame(rows)
    results.to_csv(OUT / "Stage4C_Programme_Contrasts.csv", index=False)
    # Summary of sequencing-depth difference, to aid interpretation.
    dep = focus.groupby("group").total_counts_shared_genes.agg(
        ["size", "median", "mean", "min", "max"]
    )
    dep.to_csv(OUT / "Stage4C_Library_Depth_By_Group.csv")
    log("\nLIBRARY SIZE BY PATHOLOGY GROUP")
    log(dep.round(1).to_string())
    log("\nPROGRAMME COMPARISONS (DCIS minus mixed invasive-associated regions)")
    log(
        results[
            [
                "programme",
                "DCIS_median",
                "mixed_invasive_median",
                "difference_DCIS_minus_mixed",
                "rho_with_log_depth_all_spots",
            ]
        ].to_string(index=False)
    )

    # Check direction of contrasts in each within-sample depth quartile.
    focus["depth_quartile"] = pd.qcut(
        focus.total_counts_shared_genes.rank(method="first"), 4, labels=[1, 2, 3, 4]
    ).astype(int)
    sensitivity = []
    for q, sub in focus.groupby("depth_quartile"):
        for p in PROGRAMMES:
            a = sub.loc[sub.group == "DCIS", p]
            b = sub.loc[sub.group == "Invasive-associated mixed", p]
            sensitivity.append(
                {
                    "quartile": int(q),
                    "programme": p,
                    "DCIS_n": len(a),
                    "mixed_invasive_n": len(b),
                    "difference_DCIS_minus_mixed": round(med(a) - med(b), 4)
                    if len(a) and len(b)
                    else np.nan,
                }
            )
    sens = pd.DataFrame(sensitivity)
    sens.to_csv(OUT / "Stage4C_Depth_Quartile_Sensitivity.csv", index=False)
    log("\nDEPTH-QUARTILE SENSITIVITY (DCIS minus mixed invasive-associated)")
    log(sens.to_string(index=False))

    # Spatial array map: pathology region categories (not individual cells).
    fig, ax = plt.subplots(figsize=(8.5, 7))
    specs = [("DCIS", "tab:orange"), ("Invasive cancer + lymphocytes", "tab:blue")]
    other = ~df.pathology.isin(GROUPS)
    ax.scatter(
        df.loc[other, "array_col"],
        df.loc[other, "array_row"],
        s=12,
        marker="h",
        color="lightgray",
        label="Other regions",
        rasterized=True,
    )
    for name, color in specs:
        sub = df.loc[df.pathology.eq(name)]
        ax.scatter(
            sub.array_col,
            sub.array_row,
            s=21,
            marker="h",
            color=color,
            label=f"{name} (n={len(sub)})",
            rasterized=True,
        )
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Visium array column", fontweight="bold")
    ax.set_ylabel("Visium array row", fontweight="bold")
    ax.set_title(f"{SAMPLE} | Pathology regions", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(
        FIG / "CID44971_DCIS_Invasive_Pathology_Map.webp", dpi=260, bbox_inches="tight"
    )
    plt.close(fig)

    # Signature contrast across pathology categories, no statistical claims.
    fig, ax = plt.subplots(figsize=(10, 5.5))
    xs = np.arange(len(results))
    w = 0.37
    ax.bar(xs - w / 2, results.DCIS_median, width=w, label="DCIS")
    ax.bar(
        xs + w / 2,
        results.mixed_invasive_median,
        width=w,
        label="Invasive cancer + lymphocytes",
    )
    ax.set_xticks(
        xs,
        [x.replace("_associated", "").replace("_", " ") for x in results.programme],
        rotation=25,
        ha="right",
    )
    ax.set_ylabel("Median within-specimen programme score", fontweight="bold")
    ax.set_title(
        f"{SAMPLE} | Descriptive pathology programme profiles", fontweight="bold"
    )
    ax.axhline(0, color="gray", lw=0.8)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "CID44971_Programme_Medians.webp", dpi=260, bbox_inches="tight")
    plt.close(fig)

    # Show actual values spatially without presenting raw scores as cell types.
    for p in ["Proliferation_associated", "Stromal_ECM"]:
        fig, ax = plt.subplots(figsize=(8.5, 7))
        v = df[p].to_numpy()
        vmin, vmax = np.percentile(v, [2, 98])
        if vmax <= vmin:
            vmax = vmin + 1e-6
        im = ax.scatter(
            df.array_col,
            df.array_row,
            c=v,
            cmap="coolwarm",
            vmin=vmin,
            vmax=vmax,
            s=19,
            marker="h",
            edgecolors="none",
            rasterized=True,
        )
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Visium array column", fontweight="bold")
        ax.set_ylabel("Visium array row", fontweight="bold")
        ax.set_title(f"{SAMPLE} | {p.replace('_', ' ')}", fontweight="bold")
        cb = fig.colorbar(im, ax=ax, shrink=0.8)
        cb.set_label("Within-specimen programme score", fontweight="bold")
        fig.tight_layout()
        fig.savefig(FIG / f"CID44971_{p}_Map.webp", dpi=260, bbox_inches="tight")
        plt.close(fig)

    focus.to_csv(OUT / "Stage4C_Comparison_Spots.csv", index=False)
    log("\nINTERPRETATION LIMITS")
    log(
        "The mixed invasive-associated group explicitly contains lymphocytes; this is not a controlled DCIS-versus-invasive comparison."
    )
    log(
        "No independent patient replication (n=1), spatial inference, or adjustment for region-specific cell mixtures."
    )
    log(
        "Clinical conclusions about progression or invasion are NOT supported by this contrast."
    )
    log(
        "Review histology overlay before making claims about region boundaries or programme localisation."
    )
    log("\nPASS | Saved 4 standalone WebP figures, 4 CSV files, and this text summary.")
    log(f"Outputs: {OUT}")
except Exception as exc:
    log(f"ERROR: {type(exc).__name__}: {exc}")
    log(traceback.format_exc())
finally:
    text = "\n".join(REPORT)
    (OUT / "Stage4C_Consolidated_Results.txt").write_text(text, encoding="utf-8")
    print(
        "\n" + "=" * 80 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 80
    )
    print(text)

# ============================================================================
# Neighbourhood and sequencing-depth sensitivity checks
# ============================================================================

from pathlib import Path
import traceback
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors


ROOT = PROJECT_ROOT
M = ROOT / "Results" / "Stage1D"
S = ROOT / "Results" / "Stage3B"
OUT = ROOT / "Results" / "Stage6B"
OUT.mkdir(parents=True, exist_ok=True)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
PROGRAMMES = [
    "B_plasma_associated",
    "Stromal_ECM",
    "T_cell_associated",
    "Epithelial",
    "Proliferation_associated",
    "Hypoxia_associated",
]
K_VALUES = [4, 6, 10]
report, spatial_rows, gradient_rows, contrast_rows, failures = [], [], [], [], []


def log(x=""):
    report.append(str(x))


def rho(x, y):
    with np.errstate(invalid="ignore"):
        v = float(spearmanr(x, y).statistic)
    return v


def load_sample(sample):
    meta = pd.read_csv(M / f"{sample}_spot_metadata.csv", index_col="barcode")
    scores = pd.read_csv(S / f"{sample}_spatial_programme_scores.csv")
    if "barcode" not in scores.columns:
        raise ValueError("Missing barcode in score file")
    scores["barcode"] = scores["barcode"].astype(str)
    meta.index = meta.index.astype(str)
    if meta.index.has_duplicates or scores.barcode.duplicated().any():
        raise ValueError("Duplicated barcodes")
    if len(meta) != len(scores) or set(meta.index) != set(scores.barcode):
        raise ValueError("Barcode mismatch")
    missing = set(PROGRAMMES) - set(scores.columns)
    if missing:
        raise ValueError(f"Missing score columns: {sorted(missing)}")
    df = meta.join(
        scores.set_index("barcode")[PROGRAMMES], how="left", validate="one_to_one"
    )
    required = ["array_row", "array_col", "total_counts_shared_genes"] + PROGRAMMES
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="raise")
    if not np.isfinite(df[required].to_numpy(dtype=float)).all():
        raise ValueError("Nonfinite input")
    if (df.total_counts_shared_genes <= 0).any():
        raise ValueError("Nonpositive depth")
    return df


def residualize(y, depth):
    # Diagnostic removal of linear association with log library size, not a
    # claim of unbiased normalization. Original scores remain unchanged.
    d = np.column_stack([np.ones(len(depth)), depth])
    b = np.linalg.lstsq(d, y, rcond=None)[0]
    return y - d @ b


def gradient(sample, df, anchor, programme):
    group = df["Classification"].fillna("[MISSING]").astype(str).to_numpy()
    mask = group == anchor
    if mask.sum() < 10 or (~mask).sum() < 10:
        return
    xy = df[["array_row", "array_col"]].to_numpy(dtype=float)
    dist = cKDTree(xy[mask]).query(xy[~mask], k=1)[0]
    depth = np.log1p(df.total_counts_shared_genes.to_numpy(dtype=float)[~mask])
    score = df[programme].to_numpy(dtype=float)[~mask]
    adjusted = residualize(score, depth)
    near = dist < 2
    far = dist >= 10
    gradient_rows.append(
        dict(
            sample=sample,
            anchor=anchor,
            programme=programme,
            n_anchor=int(mask.sum()),
            n_other=len(score),
            n_near=int(near.sum()),
            n_far=int(far.sum()),
            rho_raw=rho(dist, score),
            rho_depth_adjusted=rho(dist, adjusted),
            rho_depth=rho(dist, depth),
            near_minus_far_raw=float(np.median(score[near]) - np.median(score[far]))
            if near.any() and far.any()
            else np.nan,
            near_minus_far_residual=float(
                np.median(adjusted[near]) - np.median(adjusted[far])
            )
            if near.any() and far.any()
            else np.nan,
        )
    )


log("PROJECT 2 | STAGE 6B — FOCUSED ROBUSTNESS / BIOLOGICAL VALIDATION")
log("Input: existing Stage1D coordinates/annotations and Stage3B signature scores.")
log(
    "Sensitivity: spatial neighbourhood size k=4,6,10 and linear log-depth residual diagnostic."
)
log(
    "No filtering, source modification, spatial inference p-values or patient-level subtype tests."
)

for sample in SAMPLES:
    try:
        df = load_sample(sample)
        xy = df[["array_row", "array_col"]].to_numpy(dtype=float)
        if len(np.unique(xy, axis=0)) != len(xy):
            raise ValueError("Duplicate spot coordinates")
        depth = np.log1p(df.total_counts_shared_genes.to_numpy(dtype=float))
        nn = NearestNeighbors(n_neighbors=max(K_VALUES) + 1).fit(xy)
        indexes = nn.kneighbors(xy, return_distance=False)[:, 1:]
        for programme in PROGRAMMES:
            score = df[programme].to_numpy(dtype=float)
            res = residualize(score, depth)
            for k in K_VALUES:
                near_mean = score[indexes[:, :k]].mean(axis=1)
                res_mean = res[indexes[:, :k]].mean(axis=1)
                spatial_rows.append(
                    dict(
                        sample=sample,
                        programme=programme,
                        k=k,
                        neighbour_r_raw=float(np.corrcoef(score, near_mean)[0, 1]),
                        neighbour_r_residual=float(np.corrcoef(res, res_mean)[0, 1]),
                        spearman_score_depth=rho(score, depth),
                    )
                )
        if sample == "1142243F":
            gradient(sample, df, "Necrosis", "Hypoxia_associated")
        if sample in ["1160920F", "CID44971", "CID4535"]:
            for p in ["T_cell_associated", "B_plasma_associated", "Epithelial"]:
                gradient(sample, df, "Lymphocytes", p)
        if sample == "CID44971":
            cls = df.Classification.fillna("[MISSING]").astype(str).to_numpy()
            a = cls == "DCIS"
            b = cls == "Invasive cancer + lymphocytes"
            if a.sum() < 10 or b.sum() < 10:
                raise ValueError("Required DCIS/invasive regions absent")
            for p in PROGRAMMES:
                val = df[p].to_numpy(dtype=float)
                resid = residualize(val, depth)
                contrast_rows.append(
                    dict(
                        sample=sample,
                        programme=p,
                        n_DCIS=int(a.sum()),
                        n_mixed_invasive=int(b.sum()),
                        DCIS_minus_mixed_raw=float(
                            np.median(val[a]) - np.median(val[b])
                        ),
                        DCIS_minus_mixed_residual=float(
                            np.median(resid[a]) - np.median(resid[b])
                        ),
                    )
                )
        log(f"{sample}: PASS | {len(df):,} spots | 6 programmes x 3 neighbour settings")
    except Exception as e:
        failures.append(dict(sample=sample, error=f"{type(e).__name__}: {e}"))
        log(f"{sample}: ERROR {type(e).__name__}: {e}")
        log(traceback.format_exc())

sp = pd.DataFrame(spatial_rows)
gr = pd.DataFrame(gradient_rows)
ct = pd.DataFrame(contrast_rows)
if not sp.empty:
    sp.to_csv(OUT / "Spatial_Neighbourhood_Sensitivity.csv", index=False)
if not gr.empty:
    gr.to_csv(OUT / "Distance_Gradient_Depth_Sensitivity.csv", index=False)
if not ct.empty:
    ct.to_csv(OUT / "DCIS_Contrast_Depth_Sensitivity.csv", index=False)
if failures:
    pd.DataFrame(failures).to_csv(OUT / "Errors.csv", index=False)
log("\n" + "=" * 90 + "\nFINAL ANALYSIS SUMMARY\n" + "=" * 90)
if not sp.empty:
    stability = sp.groupby(["sample", "programme"], as_index=False).agg(
        min_neighbour_r_raw=("neighbour_r_raw", "min"),
        max_neighbour_r_raw=("neighbour_r_raw", "max"),
        min_neighbour_r_depth_residual=("neighbour_r_residual", "min"),
        max_neighbour_r_depth_residual=("neighbour_r_residual", "max"),
        rho_score_depth=("spearman_score_depth", "first"),
    )
    stability.to_csv(OUT / "Spatial_Stability_Summary.csv", index=False)
    log("\nSPATIAL CORRELATION ROBUSTNESS (k=4,6,10)")
    log(stability.round(3).to_string(index=False))
if not gr.empty:
    log("\nFOCUSED NECROSIS / IMMUNE DISTANCE GRADIENTS (non-anchor spots)")
    log(gr.round(3).to_string(index=False))
if not ct.empty:
    log("\nDCIS VS MIXED INVASIVE REGIONS — LOG-DEPTH SENSITIVITY")
    log(ct.round(3).to_string(index=False))
log(
    f"\nSpecimens completed: {len(SAMPLES) - len(failures)}/{len(SAMPLES)}; errors: {len(failures)}"
)
log(
    "INTERPRETATION: residuals remove only a linear log-depth association; they may remove real biological variation."
)
log(
    "Spatial neighbour correlations and anatomical gradients remain descriptive, with spatial dependence."
)
log(
    "Depth robustness alone does not establish tumour-intrinsic biology or causal spatial mechanisms."
)
log(
    "No claims about patient-level replication for single-specimen comparisons or immune infiltration into tumour cells."
)
log(f"Outputs: {OUT}")
result = "\n".join(report)
(OUT / "Stage6B_Consolidated_Report.txt").write_text(result, encoding="utf-8")
print(result)
if failures:
    raise RuntimeError(f"Spatial analysis incomplete: {len(failures)} sample(s) failed; see Stage6B_Consolidated_Report.txt")
