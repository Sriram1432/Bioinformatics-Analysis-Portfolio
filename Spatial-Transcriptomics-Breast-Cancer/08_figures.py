"""Regenerate the fourteen project figures.

Figures are assembled from saved numerical outputs, source histology and
stored spatial coordinates. The script writes high-resolution WebP and SVG
exports while preserving the established figure order and transformations.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ============================================================================
# Dataset histopathology landscape
# ============================================================================

from io import BytesIO
import json
from pathlib import Path
import tarfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


SPATIAL_ARCHIVE = PROJECT_ROOT / "Data" / "spatial.tar.gz"
STAGE1D = PROJECT_ROOT / "Results" / "Stage1D"
OUTPUT_DIR = PROJECT_ROOT / "Results" / "Scientific_Figure_Atlas"

SAMPLES = [
    "1142243F",
    "1160920F",
    "CID4290",
    "CID4465",
    "CID44971",
    "CID4535",
]

# Exact source labels are assigned stable colours across specimens. Related
# tissue descriptions share a colour family, while distinct annotations remain
# visually separable. Black marker edges provide additional discrimination on
# heterogeneous H&E backgrounds.
PATHOLOGY_COLOURS = {
    "Invasive cancer": "#0057B8",
    "Invasive cancer + lymphocytes": "#00A6A6",
    "Invasive cancer + stroma": "#6A3D9A",
    "Invasive cancer + stroma + lymphocytes": "#D81B60",
    "Invasive cancer + adipose tissue + lymphocytes": "#4ECDC4",
    "Cancer trapped in lymphocyte aggregation": "#8C564B",
    "DCIS": "#F28E2B",
    "Necrosis": "#3B2F2F",
    "Stroma": "#2CA02C",
    "Stroma + adipose tissue": "#8DAA45",
    "Lymphocytes": "#E15759",
    "TLS": "#B07AA1",
    "Normal glands + lymphocytes": "#EDC948",
    "Normal + stroma + lymphocytes": "#C49A00",
    "Normal duct": "#FFBE0B",
    "Adipose tissue": "#56B4E9",
    "Uncertain": "#7F7F7F",
    "[MISSING]": "#D9D9D9",
}

LEGEND_ORDER = list(PATHOLOGY_COLOURS)


def read_histology_and_scale(
    archive: tarfile.TarFile, sample: str
) -> tuple[np.ndarray, float]:
    """Load the supplied high-resolution H&E and Space Ranger scale factor."""
    prefix = f"spatial/{sample}_spatial/"
    image_member = archive.extractfile(prefix + "tissue_hires_image.png")
    scale_member = archive.extractfile(prefix + "scalefactors_json.json")
    if image_member is None or scale_member is None:
        raise FileNotFoundError(f"{sample}: histology image or scale factors missing")
    image = np.asarray(Image.open(BytesIO(image_member.read())).convert("RGB"))
    factors = json.loads(scale_member.read().decode("utf-8"))
    scale = float(factors["tissue_hires_scalef"])
    if image.ndim != 3 or image.shape[2] != 3 or not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"{sample}: invalid image or scale factor")
    return image, scale


def read_spots(sample: str) -> pd.DataFrame:
    """Read retained spots and validate the fields used in the figure."""
    path = STAGE1D / f"{sample}_spot_metadata.csv"
    spots = pd.read_csv(path)
    required = {"barcode", "subtype", "Classification", "pixel_row", "pixel_col"}
    missing = required.difference(spots.columns)
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
    if spots["barcode"].duplicated().any():
        raise ValueError(f"{sample}: duplicated retained-spot barcodes")
    for column in ("pixel_row", "pixel_col"):
        spots[column] = pd.to_numeric(spots[column], errors="raise")
    if not np.isfinite(spots[["pixel_row", "pixel_col"]].to_numpy(float)).all():
        raise ValueError(f"{sample}: non-finite spatial coordinates")
    subtypes = spots["subtype"].dropna().astype(str).unique()
    if len(subtypes) != 1:
        raise ValueError(f"{sample}: expected one verified subtype, found {subtypes}")
    spots["Classification"] = spots["Classification"].fillna("[MISSING]").astype(str)
    unknown = set(spots["Classification"]).difference(PATHOLOGY_COLOURS)
    if unknown:
        raise ValueError(f"{sample}: unmapped pathology labels {sorted(unknown)}")
    return spots


def crop_bounds(
    x: np.ndarray, y: np.ndarray, width: int, height: int
) -> tuple[int, int, int, int]:
    """Crop around all analysed spots while retaining undistorted image geometry."""
    padding = max(55.0, 0.055 * max(float(np.ptp(x)), float(np.ptp(y))))
    left = max(0, int(np.floor(x.min() - padding)))
    right = min(width, int(np.ceil(x.max() + padding)))
    top = max(0, int(np.floor(y.min() - padding)))
    bottom = min(height, int(np.ceil(y.max() + padding)))
    if right <= left or bottom <= top:
        raise ValueError("Invalid histology crop")
    return left, right, top, bottom


def build_figure() -> tuple[Path, Path, list[str]]:
    """Create WebP and SVG outputs and return validation notes."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "text.color": "black",
            "axes.labelcolor": "black",
            "svg.fonttype": "none",
        }
    )

    figure, axes = plt.subplots(2, 3, figsize=(20.0, 13.4), facecolor="white")
    # The six tissue panels occupy the full width. The exact-category legend is
    # placed below them in three columns to avoid a large unused right margin.
    figure.subplots_adjust(
        left=0.015, right=0.985, bottom=0.170, top=0.985, wspace=0.025, hspace=0.035
    )
    observed_labels: set[str] = set()
    checks: list[str] = []

    with tarfile.open(SPATIAL_ARCHIVE, "r:gz") as archive:
        for axis, sample in zip(axes.flat, SAMPLES):
            spots = read_spots(sample)
            image, scale = read_histology_and_scale(archive, sample)
            x = spots["pixel_col"].to_numpy(float) * scale
            y = spots["pixel_row"].to_numpy(float) * scale
            inside = (x >= 0) & (y >= 0) & (x < image.shape[1]) & (y < image.shape[0])
            if not inside.all():
                raise ValueError(
                    f"{sample}: {int((~inside).sum())} spot centres outside H&E"
                )

            left, right, top, bottom = crop_bounds(
                x, y, width=image.shape[1], height=image.shape[0]
            )
            axis.imshow(image, origin="upper", interpolation="nearest")

            labels = spots["Classification"].to_numpy(str)
            for label in LEGEND_ORDER:
                selected = labels == label
                if not selected.any():
                    continue
                observed_labels.add(label)
                axis.scatter(
                    x[selected],
                    y[selected],
                    s=6.2,
                    marker="o",
                    c=PATHOLOGY_COLOURS[label],
                    alpha=0.40,
                    linewidths=0.14,
                    edgecolors="black",
                    rasterized=True,
                )

            axis.set_xlim(left, right)
            axis.set_ylim(bottom, top)
            axis.set_aspect("equal", adjustable="box")
            axis.set_xticks([])
            axis.set_yticks([])
            for spine in axis.spines.values():
                spine.set_color("black")
                spine.set_linewidth(0.7)

            subtype = spots["subtype"].dropna().astype(str).iloc[0]
            axis.text(
                0.018,
                0.975,
                f"{sample} · {subtype}\n{len(spots):,} retained spots",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=10.5,
                fontweight="bold",
                color="black",
                bbox={
                    "boxstyle": "round,pad=0.28",
                    "facecolor": "white",
                    "edgecolor": "black",
                    "linewidth": 0.55,
                    "alpha": 0.88,
                },
            )
            checks.append(
                f"{sample}: image={image.shape[1]}x{image.shape[0]}, "
                f"scale={scale:.9g}, inside=100.0%, spots={len(spots)}"
            )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=7,
            markerfacecolor=PATHOLOGY_COLOURS[label],
            markeredgecolor="black",
            markeredgewidth=0.35,
            label=label,
        )
        for label in LEGEND_ORDER
        if label in observed_labels
    ]
    legend = figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.008),
        frameon=False,
        title="Source pathology annotation",
        title_fontproperties={"weight": "bold", "size": 11.5},
        prop={"size": 9.2},
        ncol=3,
        columnspacing=2.1,
        labelspacing=0.72,
        handletextpad=0.65,
        borderaxespad=0,
    )
    for text in legend.get_texts():
        text.set_color("black")

    webp = OUTPUT_DIR / "Figure01_Dataset_Histopathological_Landscape.webp"
    svg = OUTPUT_DIR / "Figure01_Dataset_Histopathological_Landscape.svg"
    figure.savefig(
        webp, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.08
    )
    figure.savefig(svg, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(figure)
    return webp, svg, checks


if __name__ == "__main__":
    webp_path, svg_path, validation = build_figure()
    print("Figure 1 generated from the source H&E image and Stage1D metadata.")
    for item in validation:
        print(item)
    print(webp_path)
    print(svg_path)

# ============================================================================
# Molecular, spatial and pathway compositions
# ============================================================================

from io import BytesIO
import json
from pathlib import Path
import tarfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr


ROOT = PROJECT_ROOT
RESULTS = ROOT / "Results"
OUT = RESULTS / "Scientific_Figure_Atlas"
SPATIAL_ARCHIVE = ROOT / "Data" / "spatial.tar.gz"
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
SUBTYPE = {
    "1142243F": "TNBC",
    "1160920F": "TNBC",
    "CID4290": "ER",
    "CID4465": "TNBC",
    "CID44971": "TNBC",
    "CID4535": "ER",
}
PROGRAMMES = [
    "Epithelial",
    "T_cell_associated",
    "B_plasma_associated",
    "Stromal_ECM",
    "Hypoxia_associated",
    "Proliferation_associated",
]
PROGRAMME_LABEL = {
    "Epithelial": "Epithelial",
    "T_cell_associated": "T-cell associated",
    "B_plasma_associated": "B/plasma associated",
    "Stromal_ECM": "Stromal ECM",
    "Hypoxia_associated": "Hypoxia associated",
    "Proliferation_associated": "Proliferation associated",
}
PATH_COLORS = {
    "Invasive cancer": "#0057B8",
    "Invasive cancer + lymphocytes": "#00A6A6",
    "Invasive cancer + stroma": "#6A3D9A",
    "Invasive cancer + stroma + lymphocytes": "#D81B60",
    "Invasive cancer + adipose tissue + lymphocytes": "#4ECDC4",
    "Cancer trapped in lymphocyte aggregation": "#8C564B",
    "DCIS": "#F28E2B",
    "Necrosis": "#3B2F2F",
    "Stroma": "#2CA02C",
    "Stroma + adipose tissue": "#8DAA45",
    "Lymphocytes": "#E15759",
    "TLS": "#B07AA1",
    "Normal glands + lymphocytes": "#EDC948",
    "Normal + stroma + lymphocytes": "#C49A00",
    "Normal duct": "#FFBE0B",
    "Adipose tissue": "#56B4E9",
    "Uncertain": "#7F7F7F",
    "[MISSING]": "#D9D9D9",
}
SAMPLE_COLORS = {
    s: c
    for s, c in zip(
        SAMPLES, ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#6A3D9A"]
    )
}


def setup_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "svg.fonttype": "none",
        }
    )


def save(fig, number: int, slug: str) -> tuple[Path, Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    webp = OUT / f"Figure{number:02d}_{slug}.webp"
    svg = OUT / f"Figure{number:02d}_{slug}.svg"
    # 320 dpi keeps even tightly cropped single-panel figures above the
    # atlas target of 3,000 pixels in width.
    fig.savefig(webp, dpi=320, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(svg, facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return webp, svg


def inset(ax, text, x=0.018, y=0.975, size=10):
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=size,
        fontweight="bold",
        color="black",
        bbox={
            "boxstyle": "round,pad=.25",
            "facecolor": "white",
            "edgecolor": "black",
            "linewidth": 0.5,
            "alpha": 0.86,
        },
    )


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.8)


def load_meta(sample):
    d = pd.read_csv(RESULTS / "Stage1D" / f"{sample}_spot_metadata.csv")
    req = {
        "barcode",
        "Classification",
        "array_row",
        "array_col",
        "pixel_row",
        "pixel_col",
        "subtype",
    }
    if req - set(d):
        raise ValueError(f"{sample} metadata missing {sorted(req - set(d))}")
    if d.barcode.duplicated().any():
        raise ValueError(f"{sample}: duplicate barcodes")
    d["Classification"] = d.Classification.fillna("[MISSING]").astype(str)
    return d


def load_scores(sample):
    d = pd.read_csv(RESULTS / "Stage3B" / f"{sample}_spatial_programme_scores.csv")
    if {"barcode", *PROGRAMMES} - set(d):
        raise ValueError(f"{sample}: incomplete programme scores")
    return d


def merged_meta_scores(sample):
    a, b = load_meta(sample), load_scores(sample)
    if len(a) != len(b) or set(a.barcode) != set(b.barcode):
        raise ValueError(f"{sample}: barcode mismatch")
    return a.merge(b[["barcode", *PROGRAMMES]], on="barcode", validate="one_to_one")


def histology(tf, sample):
    prefix = f"spatial/{sample}_spatial/"
    imh = tf.extractfile(prefix + "tissue_hires_image.png")
    sfh = tf.extractfile(prefix + "scalefactors_json.json")
    if imh is None or sfh is None:
        raise FileNotFoundError(sample)
    image = np.asarray(Image.open(BytesIO(imh.read())).convert("RGB"))
    scale = float(json.loads(sfh.read().decode())["tissue_hires_scalef"])
    return image, scale


def spatial_geometry(df, image, scale):
    x = df.pixel_col.to_numpy(float) * scale
    y = df.pixel_row.to_numpy(float) * scale
    inside = (x >= 0) & (y >= 0) & (x < image.shape[1]) & (y < image.shape[0])
    if not inside.all():
        raise ValueError(f"{(~inside).sum()} spots outside H&E")
    pad = max(55, 0.055 * max(np.ptp(x), np.ptp(y)))
    bounds = (
        max(0, int(x.min() - pad)),
        min(image.shape[1], int(x.max() + pad)),
        max(0, int(y.min() - pad)),
        min(image.shape[0], int(y.max() + pad)),
    )
    return x, y, bounds


def show_hist(ax, image, bounds):
    left, right, top, bottom = bounds
    ax.imshow(image, origin="upper", interpolation="nearest")
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.7)
        s.set_color("black")


def score_overlay(fig, ax, image, x, y, bounds, values, label, rings=None):
    show_hist(ax, image, bounds)
    lim = float(np.percentile(np.abs(values), 98))
    lim = max(lim, 1e-6)
    art = ax.scatter(
        x,
        y,
        c=values,
        s=7,
        alpha=0.54,
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
        linewidths=0,
        rasterized=True,
    )
    if rings:
        labels = rings["labels"]
        for name, color in rings["targets"]:
            m = labels == name
            ax.scatter(
                x[m],
                y[m],
                s=17,
                facecolors="none",
                edgecolors=color,
                linewidths=0.55,
                alpha=0.9,
                rasterized=True,
            )
    inset(ax, label)
    cb = fig.colorbar(art, ax=ax, shrink=0.70, pad=0.012)
    cb.set_label("Within-tumour mean gene z-score", fontweight="bold", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)


def heatmap(ax, matrix, xlabels, ylabels, cmap, norm, annotate=True):
    im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(xlabels)), xlabels, rotation=38, ha="right")
    ax.set_yticks(range(len(ylabels)), ylabels)
    if annotate:
        threshold = (np.nanmin(matrix) + np.nanmax(matrix)) / 2
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                v = matrix[i, j]
                if np.isfinite(v):
                    ax.text(
                        j,
                        i,
                        f"{v:.2f}",
                        ha="center",
                        va="center",
                        fontsize=7.5,
                        color="white" if v > threshold else "black",
                        fontweight="bold",
                    )
    return im


def figure02():
    chosen = ["1160920F", "CID44971"]
    fig = plt.figure(figsize=(18, 9.7), facecolor="white")
    grid = fig.add_gridspec(
        2,
        4,
        width_ratios=[1, 1, 1, 0.72],
        left=0.055,
        right=0.985,
        top=0.98,
        bottom=0.09,
        wspace=0.26,
        hspace=0.28,
    )
    axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(3)] for i in range(2)])
    legend_axes = [fig.add_subplot(grid[i, 3]) for i in range(2)]
    cluster_palette = [
        "#0072B2",
        "#E69F00",
        "#009E73",
        "#CC79A7",
        "#D55E00",
        "#56B4E9",
        "#F0E442",
        "#6A3D9A",
        "#7F7F7F",
        "#8C564B",
    ]
    for row, sample in enumerate(chosen):
        meta = load_meta(sample).set_index("barcode")
        pca = pd.read_csv(
            RESULTS / "Stage2B" / f"{sample}_PCA30_spot_scores.csv", index_col="barcode"
        )
        emb = pd.read_csv(
            RESULTS / "Stage2E" / f"{sample}_log_zscore_spot_clusters_UMAP.csv"
        ).set_index("barcode")
        if not meta.index.equals(pca.index) or not meta.index.equals(emb.index):
            raise ValueError("Figure02 order mismatch")
        path = meta.Classification.to_numpy(str)
        ax = axes[row, 0]
        for label in pd.unique(path):
            m = path == label
            ax.scatter(
                pca.PC1[m],
                pca.PC2[m],
                s=5,
                alpha=0.52,
                c=PATH_COLORS[label],
                linewidths=0,
                rasterized=True,
            )
        ax.set_xlabel("PC1", fontweight="bold")
        ax.set_ylabel("PC2", fontweight="bold")
        inset(ax, f"{sample} · {SUBTYPE[sample]}\nPCA · pathology")
        clean(ax)
        ax = axes[row, 1]
        clusters = np.sort(emb.cluster.unique().astype(int))
        if not np.array_equal(clusters, np.arange(len(clusters))):
            raise ValueError(f"{sample}: nonconsecutive clusters")
        cmap = ListedColormap(cluster_palette[: len(clusters)])
        norm = BoundaryNorm(np.arange(-0.5, len(clusters) + 0.5), cmap.N)
        art = ax.scatter(
            emb.UMAP1,
            emb.UMAP2,
            c=emb.cluster,
            s=5,
            alpha=0.72,
            cmap=cmap,
            norm=norm,
            linewidths=0,
            rasterized=True,
        )
        ax.set_xlabel("UMAP 1 (embedding)", fontweight="bold")
        ax.set_ylabel("UMAP 2 (embedding)", fontweight="bold")
        inset(ax, "UMAP · transcriptional clusters")
        clean(ax)
        cb = fig.colorbar(art, ax=ax, shrink=0.58, pad=0.01, ticks=clusters)
        cb.set_label("Cluster", fontweight="bold")
        cb.ax.tick_params(labelsize=8)
        ax = axes[row, 2]
        for label in pd.unique(path):
            m = path == label
            ax.scatter(
                emb.UMAP1[m],
                emb.UMAP2[m],
                s=5,
                alpha=0.58,
                c=PATH_COLORS[label],
                linewidths=0,
                rasterized=True,
            )
        ax.set_xlabel("UMAP 1 (embedding)", fontweight="bold")
        ax.set_ylabel("UMAP 2 (embedding)", fontweight="bold")
        inset(ax, "UMAP · pathology")
        clean(ax)
        labels = [x for x in PATH_COLORS if x in set(path)]
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markersize=6,
                markerfacecolor=PATH_COLORS[x],
                markeredgecolor="black",
                markeredgewidth=0.3,
                label=x,
            )
            for x in labels
        ]
        lax = legend_axes[row]
        lax.axis("off")
        lax.legend(
            handles=handles,
            loc="center left",
            frameon=False,
            ncol=1,
            title=f"{sample}\nsource pathology",
            title_fontproperties={"weight": "bold", "size": 10},
            prop={"size": 8.5},
            handletextpad=0.55,
            labelspacing=0.7,
            borderaxespad=0,
        )
    return save(fig, 2, "Transcriptional_Heterogeneity_Pathology")


def figure03():
    d = pd.read_csv(RESULTS / "Stage3B" / "Stage3B_cluster_programme_scores.csv")
    req = {"sample", "programme", "cluster", "n_spots", "mean_score"}
    if req - set(d):
        raise ValueError("Figure03 columns")
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), facecolor="white")
    fig.subplots_adjust(
        left=0.06, right=0.98, top=0.97, bottom=0.12, wspace=0.34, hspace=0.43
    )
    for ax, sample in zip(axes.flat, SAMPLES):
        q = (
            d[d["sample"] == sample]
            .pivot(index="cluster", columns="programme", values="mean_score")
            .reindex(columns=PROGRAMMES)
        )
        a = q.to_numpy(float)
        lim = max(float(np.nanpercentile(np.abs(a), 98)), 0.05)
        im = heatmap(
            ax,
            a,
            [PROGRAMME_LABEL[x] for x in PROGRAMMES],
            [f"Cluster {int(x)}" for x in q.index],
            "coolwarm",
            TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
            True,
        )
        ax.set_xlabel("Curated expression programme", fontweight="bold")
        ax.set_ylabel(f"{sample} · {SUBTYPE[sample]}", fontweight="bold")
        cb = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.015)
        cb.set_label("Mean score\n(within tumour)", fontweight="bold", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    return save(fig, 3, "Molecular_Programmes_Transcriptional_Regions")


def figure04():
    d = pd.read_csv(RESULTS / "Stage3B" / "Stage3B_programme_summary.csv")
    q = d.pivot(
        index="sample", columns="programme", values="spatial_neighbour_r"
    ).reindex(index=SAMPLES, columns=PROGRAMMES)
    fig, ax = plt.subplots(figsize=(12.5, 7.2), facecolor="white")
    fig.subplots_adjust(left=0.18, right=0.91, bottom=0.23, top=0.97)
    im = heatmap(
        ax,
        q.to_numpy(float),
        [PROGRAMME_LABEL[x] for x in PROGRAMMES],
        [f"{s} · {SUBTYPE[s]}" for s in SAMPLES],
        "viridis",
        Normalize(0, 1),
        True,
    )
    ax.set_xlabel("Curated expression programme", fontweight="bold")
    ax.set_ylabel("Tumour specimen", fontweight="bold")
    cb = fig.colorbar(im, ax=ax, pad=0.025)
    cb.set_label("Spatial-neighbour correlation (r)", fontweight="bold")
    return save(fig, 4, "Cross_Tumour_Spatial_Organisation")


def figure05():
    sample = "1142243F"
    df = merged_meta_scores(sample)
    labels = df.Classification.to_numpy(str)
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.4), facecolor="white")
    fig.subplots_adjust(left=0.02, right=0.96, bottom=0.03, top=0.98, wspace=0.08)
    with tarfile.open(SPATIAL_ARCHIVE, "r:gz") as tf:
        image, scale = histology(tf, sample)
    x, y, bounds = spatial_geometry(df, image, scale)
    show_hist(axes[0], image, bounds)
    axes[0].scatter(x, y, s=5, c="#777777", alpha=0.12, linewidths=0, rasterized=True)
    nec = labels == "Necrosis"
    axes[0].scatter(
        x[nec],
        y[nec],
        s=11,
        c=PATH_COLORS["Necrosis"],
        alpha=0.72,
        linewidths=0.15,
        edgecolors="black",
        rasterized=True,
    )
    inset(axes[0], f"{sample} · pathology\nNecrosis (n={nec.sum():,})")
    score_overlay(
        fig,
        axes[1],
        image,
        x,
        y,
        bounds,
        df.Hypoxia_associated.to_numpy(float),
        "Hypoxia-associated score",
        {"labels": labels, "targets": [("Necrosis", PATH_COLORS["Necrosis"])]},
    )
    return save(fig, 5, "Hypoxia_Necrosis_Histology")


def figure06():
    d = pd.read_csv(RESULTS / "Stage4A" / "1142243F_signed_distance_and_scores.csv")
    req = {"is_necrosis", "signed_necrosis_distance_array", "Hypoxia_associated"}
    if req - set(d):
        raise ValueError("Figure06 columns")
    non = d[~d.is_necrosis.astype(bool)].copy()
    x = non.signed_necrosis_distance_array.to_numpy(float)
    y = non.Hypoxia_associated.to_numpy(float)
    rho = float(spearmanr(x, y).statistic)
    if not np.isclose(rho, -0.3481, atol=0.002):
        raise ValueError(f"Unexpected hypoxia rho {rho}")
    edges = [0, 2, 5, 10, np.inf]
    labels = ["<2", "2–<5", "5–<10", "≥10"]
    bins = pd.cut(x, edges, right=False, include_lowest=True, labels=labels)
    summary = (
        pd.DataFrame({"x": x, "y": y, "bin": bins})
        .groupby("bin", observed=False)
        .agg(x=("x", "median"), y=("y", "median"), n=("y", "size"))
    )
    fig, ax = plt.subplots(figsize=(12.5, 7.5), facecolor="white")
    ax.scatter(
        x,
        y,
        s=8,
        c="#6C757D",
        alpha=0.16,
        linewidths=0,
        rasterized=True,
        label="Non-necrotic Visium spots",
    )
    ax.plot(
        summary.x,
        summary.y,
        color="#D55E00",
        linewidth=2.4,
        marker="o",
        markersize=7,
        label="Distance-band median",
    )
    for r in summary.itertuples():
        ax.annotate(
            f"n={r.n:,}",
            (r.x, r.y),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )
    ax.axhline(0, color="black", linewidth=0.7, alpha=0.55)
    ax.set_xlabel(
        "Distance to nearest necrosis-annotated spot centre (Visium array units)",
        fontweight="bold",
    )
    ax.set_ylabel("Hypoxia-associated score (within tumour)", fontweight="bold")
    inset(
        ax, f"1142243F · non-necrotic spots\nSpearman ρ = {rho:.3f} · n = {len(non):,}"
    )
    ax.legend(frameon=False, loc="upper right")
    clean(ax)
    return save(fig, 6, "Hypoxia_Necrosis_Distance_Gradient")


def figure07():
    samples = ["1160920F", "CID44971", "CID4535"]
    progs = ["T_cell_associated", "B_plasma_associated"]
    fig, axes = plt.subplots(3, 2, figsize=(15.5, 18), facecolor="white")
    fig.subplots_adjust(
        left=0.02, right=0.96, bottom=0.025, top=0.985, wspace=0.10, hspace=0.05
    )
    with tarfile.open(SPATIAL_ARCHIVE, "r:gz") as tf:
        for i, sample in enumerate(samples):
            df = merged_meta_scores(sample)
            image, scale = histology(tf, sample)
            x, y, bounds = spatial_geometry(df, image, scale)
            labels = df.Classification.to_numpy(str)
            rings = {
                "labels": labels,
                "targets": [("Lymphocytes", PATH_COLORS["Lymphocytes"])],
            }
            for j, p in enumerate(progs):
                score_overlay(
                    fig,
                    axes[i, j],
                    image,
                    x,
                    y,
                    bounds,
                    df[p].to_numpy(float),
                    f"{sample} · {PROGRAMME_LABEL[p]}\nLymphocytes n={(labels == 'Lymphocytes').sum():,}",
                    rings,
                )
    return save(fig, 7, "Immune_Associated_Molecular_Architecture")


def figure08():
    bands = pd.read_csv(RESULTS / "Stage4B" / "Stage4B_immune_distance_bands.csv")
    robust = pd.read_csv(
        RESULTS / "Stage6B" / "Distance_Gradient_Depth_Sensitivity.csv"
    )
    samples = ["1160920F", "CID44971", "CID4535"]
    progs = ["T_cell_associated", "B_plasma_associated"]
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.7), facecolor="white")
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.18, top=0.97, wspace=0.25)
    x = np.arange(4)
    xt = ["<2", "2–<5", "5–<10", "≥10"]
    offsets = {"1160920F": -0.18, "CID44971": 0, "CID4535": 0.18}
    markers = {"1160920F": "o", "CID44971": "s", "CID4535": "D"}
    for ax, p in zip(axes, progs):
        for sample in samples:
            g = bands[(bands["sample"] == sample) & (bands.programme == p)]
            if len(g) != 4:
                raise ValueError(f"Figure08 {sample} {p}")
            r = robust[(robust["sample"] == sample) & (robust.programme == p)].iloc[0]
            ax.scatter(
                x + offsets[sample],
                g.median_score.to_numpy(float),
                s=74,
                marker=markers[sample],
                color=SAMPLE_COLORS[sample],
                edgecolors="black",
                linewidths=0.45,
                label=f"{sample}  spot-level ρ={r.rho_raw:+.3f}",
                zorder=3,
            )
        ax.axhline(0, color="black", linewidth=0.7, alpha=0.55)
        ax.set_xticks(x, xt)
        ax.set_xlabel(
            "Distance to lymphocyte-annotated spot (array units)", fontweight="bold"
        )
        ax.set_ylabel("Distance-bin median within-tumour score", fontweight="bold")
        inset(ax, PROGRAMME_LABEL[p])
        ax.legend(
            frameon=False,
            fontsize=8.6,
            title="Spearman ρ uses individual non-anchor spots",
            title_fontproperties={"size": 8.5, "weight": "bold"},
        )
        clean(ax)
    return save(fig, 8, "Immune_Associated_Distance_Gradients")


def figure09():
    sample = "CID44971"
    df = merged_meta_scores(sample)
    labels = df.Classification.to_numpy(str)
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.4), facecolor="white")
    fig.subplots_adjust(left=0.02, right=0.96, bottom=0.10, top=0.98, wspace=0.09)
    with tarfile.open(SPATIAL_ARCHIVE, "r:gz") as tf:
        image, scale = histology(tf, sample)
    x, y, bounds = spatial_geometry(df, image, scale)
    targets = [
        ("DCIS", PATH_COLORS["DCIS"]),
        ("Invasive cancer + lymphocytes", PATH_COLORS["Invasive cancer + lymphocytes"]),
    ]
    rings = {"labels": labels, "targets": targets}
    for ax, p in zip(axes, ["Proliferation_associated", "Stromal_ECM"]):
        score_overlay(
            fig,
            ax,
            image,
            x,
            y,
            bounds,
            df[p].to_numpy(float),
            f"{sample} · {PROGRAMME_LABEL[p]}",
            rings,
        )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor=c,
            markeredgewidth=1.5,
            label=f"{n} (n={(labels == n).sum():,})",
        )
        for n, c in targets
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2,
        frameon=False,
        prop={"size": 9},
    )
    return save(fig, 9, "DCIS_Associated_Molecular_Heterogeneity")


def figure10():
    d = pd.read_csv(RESULTS / "Stage6B" / "DCIS_Contrast_Depth_Sensitivity.csv")
    d = d.set_index("programme").reindex(PROGRAMMES).reset_index()
    if d[["DCIS_minus_mixed_raw", "DCIS_minus_mixed_residual"]].isna().any().any():
        raise ValueError("Figure10 missing")
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(12.8, 7.2), facecolor="white")
    for i, r in d.iterrows():
        ax.plot(
            [r.DCIS_minus_mixed_raw, r.DCIS_minus_mixed_residual],
            [i, i],
            color="#A7A7A7",
            linewidth=1.5,
            zorder=1,
        )
    ax.scatter(d.DCIS_minus_mixed_raw, y, s=78, c="#0072B2", label="Observed", zorder=3)
    ax.scatter(
        d.DCIS_minus_mixed_residual,
        y,
        s=78,
        c="#D55E00",
        marker="D",
        label="After linear log-depth residualisation",
        zorder=3,
    )
    ax.axvline(0, color="black", linewidth=0.9)
    ax.set_yticks(y, [PROGRAMME_LABEL[x] for x in d.programme])
    ax.invert_yaxis()
    ax.set_xlabel(
        "Median-score contrast: DCIS − Invasive cancer + lymphocytes", fontweight="bold"
    )
    ax.set_ylabel("Curated expression programme", fontweight="bold")
    inset(
        ax,
        "CID44971 · both estimates descriptive\nDCIS n=273 · mixed invasive-associated n=317",
    )
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.985, 0.93))
    ax.grid(axis="x", alpha=0.16)
    clean(ax)
    return save(fig, 10, "Histopathological_Programme_Contrasts")


def select_recurrent_genes():
    all_rows = []
    for s in SAMPLES:
        d = pd.read_csv(RESULTS / "Stage5A" / f"{s}_MoranI_gene_results.csv")
        d["sample"] = s
        all_rows.append(d)
    d = pd.concat(all_rows, ignore_index=True)
    elig = (
        d.groupby("gene_symbol")
        .agg(
            specimens=(
                "sample",
                lambda x: x[d.loc[x.index, "FDR_lt_0_05"].astype(bool)].nunique(),
            ),
            max_abs_depth_r=("pearson_r_log_depth", lambda x: float(np.max(np.abs(x)))),
            median_I=("moran_I", "median"),
        )
        .reset_index()
    )
    elig = elig[
        (elig.specimens == 6)
        & (elig.max_abs_depth_r < 0.4)
        & ~elig.gene_symbol.str.startswith(("MT-", "RPL", "RPS"))
    ]
    priority = [
        "IGHG3",
        "IGLC3",
        "IGKC",
        "JCHAIN",
        "COL1A1",
        "COL3A1",
        "COL6A2",
        "BGN",
        "AEBP1",
        "IGFBP7",
        "TIMP1",
        "SPARC",
        "DCN",
        "LUM",
        "MYL9",
        "TAGLN",
        "APOE",
        "APOC1",
        "C1QA",
        "C1QB",
        "TYROBP",
        "FCER1G",
        "HLA-DPA1",
        "HLA-DPB1",
        "CD74",
        "IFI30",
        "CTSB",
    ]
    selected = [g for g in priority if g in set(elig.gene_symbol)]
    for g in elig.sort_values("median_I", ascending=False).gene_symbol:
        if g not in selected:
            selected.append(g)
        if len(selected) >= 25:
            break
    selected = selected[:25]
    if len(selected) < 20:
        raise ValueError(f"Only {len(selected)} recurrent genes")
    matrix = (
        d[d.gene_symbol.isin(selected)]
        .pivot(index="gene_symbol", columns="sample", values="moran_I")
        .reindex(index=selected, columns=SAMPLES)
    )
    if matrix.isna().any().any():
        raise ValueError("Figure11 incomplete matrix")
    audit = elig.set_index("gene_symbol").loc[selected].reset_index()
    audit.insert(0, "selection_rank", range(1, len(audit) + 1))
    audit.to_csv(OUT / "Figure11_Selected_Recurrent_Genes.csv", index=False)
    return matrix


def figure11():
    q = select_recurrent_genes()
    a = q.to_numpy(float)
    fig, ax = plt.subplots(figsize=(12.5, 12.2), facecolor="white")
    fig.subplots_adjust(left=0.18, right=0.91, bottom=0.12, top=0.98)
    im = heatmap(
        ax,
        a,
        SAMPLES,
        q.index.tolist(),
        "magma",
        Normalize(0, float(np.nanmax(a))),
        True,
    )
    ax.set_xlabel("Tumour specimen", fontweight="bold")
    ax.set_ylabel("Recurrent spatial gene", fontweight="bold")
    cb = fig.colorbar(im, ax=ax, pad=0.025)
    cb.set_label("Moran's I", fontweight="bold")
    return save(fig, 11, "Recurrent_Spatially_Structured_Genes")


PATHWAY_SELECTION = [
    (
        "GO_Biological_Process_2025",
        "Antigen Processing and Presentation of Exogenous Peptide Antigen via MHC Class II (GO:0019886)",
        "GO · MHC-II antigen presentation",
    ),
    (
        "GO_Biological_Process_2025",
        "Extracellular Matrix Organization (GO:0030198)",
        "GO · extracellular matrix organisation",
    ),
    (
        "GO_Biological_Process_2025",
        "Positive Regulation of Immune System Process (GO:0002684)",
        "GO · positive immune-system regulation",
    ),
    (
        "GO_Biological_Process_2025",
        "Collagen Fibril Organization (GO:0030199)",
        "GO · collagen fibril organisation",
    ),
    (
        "GO_Biological_Process_2025",
        "Cytoplasmic Translation (GO:0002181)",
        "GO · cytoplasmic translation",
    ),
    (
        "GO_Biological_Process_2025",
        "Antimicrobial Humoral Response (GO:0019730)",
        "GO · antimicrobial humoral response",
    ),
    (
        "Reactome_2022",
        "Binding And Uptake Of Ligands By Scavenger Receptors R-HSA-2173782",
        "Reactome · scavenger-receptor ligand uptake",
    ),
    (
        "Reactome_2022",
        "Interferon Gamma Signaling R-HSA-877300",
        "Reactome · interferon-γ signalling",
    ),
    (
        "Reactome_2022",
        "Extracellular Matrix Organization R-HSA-1474244",
        "Reactome · extracellular matrix organisation",
    ),
    (
        "Reactome_2022",
        "Peptide Chain Elongation R-HSA-156902",
        "Reactome · peptide-chain elongation",
    ),
    (
        "Reactome_2022",
        "Initial Triggering Of Complement R-HSA-166663",
        "Reactome · initial complement triggering",
    ),
    ("Reactome_2022", "Keratinization R-HSA-6805567", "Reactome · keratinisation"),
    (
        "Reactome_2022",
        "Assembly Of Collagen Fibrils And Other Multimeric Structures R-HSA-2022090",
        "Reactome · collagen-fibril assembly",
    ),
]


def figure12():
    rows = []
    for sample in SAMPLES:
        for library, pathway, label in PATHWAY_SELECTION:
            short = "GO_BP" if library.startswith("GO_") else "Reactome"
            d = pd.read_csv(
                RESULTS / "Stage5B" / f"{sample}_{short}_enrichment_all_tested.csv"
            )
            hit = d[d.pathway == pathway]
            if hit.empty:
                rows.append(
                    {
                        "sample": sample,
                        "library": library,
                        "pathway": pathway,
                        "label": label,
                        "tested": False,
                        "overlap": np.nan,
                        "fold": np.nan,
                        "q": np.nan,
                        "significant": False,
                    }
                )
            else:
                r = hit.iloc[0]
                rows.append(
                    {
                        "sample": sample,
                        "library": library,
                        "pathway": pathway,
                        "label": label,
                        "tested": True,
                        "overlap": int(r.overlap_count),
                        "fold": float(r.fold_enrichment),
                        "q": float(r.BH_FDR),
                        "significant": bool(r.significant),
                    }
                )
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "Figure12_Selected_Pathway_Results.csv", index=False)
    fig, ax = plt.subplots(figsize=(15.5, 9.5), facecolor="white")
    fig.subplots_adjust(left=0.31, right=0.91, bottom=0.13, top=0.98)
    labels = [x[2] for x in PATHWAY_SELECTION]
    finite = audit.loc[audit.tested & (audit.overlap > 0), "fold"]
    vmax = max(3, float(np.nanpercentile(np.log2(finite), 98)))

    def enrichment_size(fold):
        return 34 + 42 * min(vmax, max(0, float(np.log2(fold))))

    for yi, label in enumerate(labels):
        for xi, sample in enumerate(SAMPLES):
            r = audit[(audit.label == label) & (audit["sample"] == sample)].iloc[0]
            if not r.tested:
                ax.scatter(xi, yi, marker="x", s=45, c="#8A8A8A", linewidths=1.3)
            elif r.overlap == 0:
                ax.scatter(
                    xi,
                    yi,
                    marker="s",
                    s=36,
                    facecolors="white",
                    edgecolors="#8A8A8A",
                    linewidths=0.9,
                )
            else:
                ax.scatter(
                    xi,
                    yi,
                    s=enrichment_size(r.fold),
                    c="#00897B" if r.significant else "#D7DEDC",
                    alpha=1,
                    edgecolors="black" if r.significant else "#8C9693",
                    linewidths=1.0 if r.significant else 0.65,
                )
    ax.set_xticks(range(len(SAMPLES)), SAMPLES, rotation=25, ha="right")
    ax.set_yticks(range(len(labels)), labels)
    ax.invert_yaxis()
    ax.set_xlabel("Tumour specimen", fontweight="bold")
    ax.set_ylabel("Gene-set overrepresentation test", fontweight="bold")
    ax.set_xlim(-0.6, len(SAMPLES) - 0.4)
    ax.grid(color="#E8E8E8", linewidth=0.6)
    ax.set_axisbelow(True)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=8,
            markerfacecolor="#00897B",
            markeredgecolor="black",
            label="BH q<0.05",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=7,
            markerfacecolor="#D7DEDC",
            markeredgecolor="#8C9693",
            label="Tested, not significant",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="none",
            markersize=7,
            markerfacecolor="white",
            markeredgecolor="#777",
            label="Tested, zero overlap",
        ),
        Line2D(
            [0],
            [0],
            marker="x",
            linestyle="none",
            markersize=7,
            color="#777",
            label="Not eligible/tested",
        ),
    ]
    for value in [1, 2, 3]:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markersize=np.sqrt(enrichment_size(2**value)),
                markerfacecolor="white",
                markeredgecolor="black",
                label=f"log₂ fold = {value}",
            )
        )
    ax.legend(
        handles=handles,
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=4,
        fontsize=7.8,
        title="Test status and enrichment magnitude (marker area)",
        title_fontproperties={"size": 8.5, "weight": "bold"},
    )
    return save(fig, 12, "Spatial_Gene_Set_Enrichment")


def main():
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    makers = [
        figure02,
        figure03,
        figure04,
        figure05,
        figure06,
        figure07,
        figure08,
        figure09,
        figure10,
        figure11,
        figure12,
    ]
    outputs = []
    for fn in makers:
        result = fn()
        outputs.extend(result)
        print(fn.__name__, *result, sep="\n  ")
    print(f"Completed {len(makers)} figures; {len(outputs)} files.")


if __name__ == "__main__":
    main()

# Convert the twelve-composition numbering to the fourteen-figure atlas.
_atlas = PROJECT_ROOT / "Results" / "Scientific_Figure_Atlas"
_renumber = {
    12: (14, "Spatial_Gene_Set_Enrichment"),
    11: (13, "Recurrent_Spatially_Structured_Genes"),
    10: (12, "Histopathological_Programme_Contrasts"),
    9: (11, "DCIS_Associated_Molecular_Heterogeneity"),
    8: (10, "Immune_Associated_Distance_Gradients"),
    7: (9, "Immune_Associated_Molecular_Architecture"),
    6: (8, "Hypoxia_Necrosis_Distance_Gradient"),
    5: (7, "Hypoxia_Necrosis_Histology"),
    4: (6, "Cross_Tumour_Spatial_Organisation"),
    3: (5, "Molecular_Programmes_Transcriptional_Regions"),
}
for _old, (_new, _name) in _renumber.items():
    for _extension in ("webp", "svg"):
        _matches = sorted(_atlas.glob(f"Figure{_old:02d}_*.{_extension}"))
        if len(_matches) != 1:
            raise RuntimeError(
                f"Expected one Figure{_old:02d} {_extension} export, found {len(_matches)}"
            )
        _destination = _atlas / f"Figure{_new:02d}_{_name}.{_extension}"
        if _destination.exists():
            _destination.unlink()
        _matches[0].replace(_destination)

# ============================================================================
# Final transcriptional and cluster-marker compositions
# ============================================================================

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = PROJECT_ROOT
OUT = ROOT / "Results" / "Scientific_Figure_Atlas"
STAGE1D = ROOT / "Results" / "Stage1D"
STAGE2E = ROOT / "Results" / "Stage2E"
STAGE3A = ROOT / "Results" / "Stage3A"
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
CLUSTER_COLOURS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#7B61A8",
    "#E69F00",
    "#56B4E9",
    "#CC79A7",
    "#8C564B",
    "#6B6B6B",
    "#829B2F",
    "#1B9E77",
    "#E7298A",
    "#7570B3",
    "#66A61E",
    "#A6761D",
]


def setup_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.labelweight": "bold",
            "axes.edgecolor": "black",
            "axes.linewidth": 0.8,
            "xtick.color": "black",
            "ytick.color": "black",
            "text.color": "black",
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def load_stage2e(sample):
    emb = pd.read_csv(STAGE2E / f"{sample}_log_zscore_spot_clusters_UMAP.csv")
    meta = pd.read_csv(STAGE1D / f"{sample}_spot_metadata.csv")
    required_emb = {"barcode", "cluster", "UMAP1", "UMAP2"}
    required_meta = {"barcode", "array_row", "array_col"}
    if not required_emb.issubset(emb.columns) or not required_meta.issubset(
        meta.columns
    ):
        raise ValueError(f"{sample}: missing Stage 2E or Stage 1D columns")
    if emb.barcode.duplicated().any() or meta.barcode.duplicated().any():
        raise ValueError(f"{sample}: duplicate barcodes")
    data = emb.merge(
        meta[["barcode", "array_row", "array_col"]],
        on="barcode",
        how="left",
        validate="one_to_one",
    )
    if data[["array_row", "array_col", "UMAP1", "UMAP2", "cluster"]].isna().any().any():
        raise ValueError(f"{sample}: incomplete coordinates or cluster labels")
    data["cluster"] = data["cluster"].astype(int)
    return data


def specimen_label(ax, sample, n):
    ax.text(
        0.02,
        0.98,
        f"{sample}  |  n={n:,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.22",
            facecolor="white",
            edgecolor="black",
            linewidth=0.55,
            alpha=0.92,
        ),
        zorder=10,
    )


def cluster_legend(ax, data):
    counts = data.groupby("cluster").size().sort_index()
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=4.3,
            markerfacecolor=CLUSTER_COLOURS[int(cluster) % len(CLUSTER_COLOURS)],
            markeredgecolor="none",
            label=f"{int(cluster)} ({int(count):,})",
        )
        for cluster, count in counts.items()
    ]
    ax.legend(
        handles=handles,
        title="Cluster (spots)",
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=min(5, len(handles)),
        fontsize=7.1,
        title_fontproperties={"size": 7.5, "weight": "bold"},
        columnspacing=0.9,
        handletextpad=0.3,
        borderaxespad=0,
    )


def save(fig, basename):
    webp = OUT / f"{basename}.webp"
    svg = OUT / f"{basename}.svg"
    fig.savefig(webp, dpi=320, bbox_inches="tight", pad_inches=0.10, format="webp")
    fig.savefig(svg, bbox_inches="tight", pad_inches=0.10, format="svg")
    plt.close(fig)
    return webp, svg


def figure02_umap():
    fig, axes = plt.subplots(2, 3, figsize=(17.4, 11.2), facecolor="white")
    for ax, sample in zip(axes.flat, SAMPLES):
        data = load_stage2e(sample)
        colours = [CLUSTER_COLOURS[c % len(CLUSTER_COLOURS)] for c in data.cluster]
        ax.scatter(
            data.UMAP1,
            data.UMAP2,
            c=colours,
            s=5.0,
            alpha=0.80,
            linewidths=0,
            rasterized=False,
        )
        specimen_label(ax, sample, len(data))
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.grid(False)
        ax.spines[["top", "right"]].set_visible(False)
        cluster_legend(ax, data)
    fig.subplots_adjust(
        left=0.055, right=0.985, top=0.975, bottom=0.095, wspace=0.22, hspace=0.48
    )
    return save(fig, "Figure02_Transcriptional_UMAP_Atlas")


def figure03_spatial_clusters():
    fig, axes = plt.subplots(2, 3, figsize=(17.4, 11.4), facecolor="white")
    for ax, sample in zip(axes.flat, SAMPLES):
        data = load_stage2e(sample)
        colours = [CLUSTER_COLOURS[c % len(CLUSTER_COLOURS)] for c in data.cluster]
        ax.scatter(
            data.array_col,
            data.array_row,
            c=colours,
            s=7.2,
            alpha=0.92,
            linewidths=0,
            rasterized=False,
        )
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")
        specimen_label(ax, sample, len(data))
        ax.set_xlabel("Visium array column")
        ax.set_ylabel("Visium array row")
        ax.grid(False)
        ax.spines[["top", "right"]].set_visible(False)
        cluster_legend(ax, data)
    fig.subplots_adjust(
        left=0.055, right=0.985, top=0.975, bottom=0.09, wspace=0.22, hspace=0.50
    )
    return save(fig, "Figure03_Spatial_Transcriptional_Clusters")


def load_marker_matrix(sample):
    top = pd.read_csv(STAGE3A / f"{sample}_top_descriptive_markers.csv")
    full = pd.read_csv(STAGE3A / f"{sample}_all_gene_descriptive_contrasts.csv.gz")
    required_top = {"cluster", "rank", "gene_symbol"}
    required_full = {"cluster", "gene_symbol", "mean_log_cluster"}
    if not required_top.issubset(top.columns) or not required_full.issubset(
        full.columns
    ):
        raise ValueError(f"{sample}: incomplete Stage 3A marker tables")
    selected = []
    for cluster in sorted(top.cluster.unique()):
        genes = (
            top.loc[(top.cluster == cluster) & (top["rank"] <= 2)]
            .sort_values("rank")
            .gene_symbol
        )
        for gene in genes:
            if gene not in selected:
                selected.append(gene)
    matrix = full[full.gene_symbol.isin(selected)].pivot(
        index="cluster", columns="gene_symbol", values="mean_log_cluster"
    )
    matrix = matrix.reindex(index=sorted(full.cluster.unique()), columns=selected)
    if matrix.empty or matrix.isna().any().any():
        raise ValueError(f"{sample}: incomplete marker-expression matrix")
    return matrix


def figure04_marker_expression():
    fig, axes = plt.subplots(2, 3, figsize=(18.2, 12.4), facecolor="white")
    summary = pd.read_csv(STAGE3A / "Stage3A_sample_summary.csv").set_index("sample")
    for ax, sample in zip(axes.flat, SAMPLES):
        matrix = load_marker_matrix(sample)
        image = ax.imshow(
            matrix.to_numpy(float),
            aspect="auto",
            interpolation="nearest",
            cmap="viridis",
        )
        ax.set_xticks(
            np.arange(matrix.shape[1]),
            matrix.columns,
            rotation=58,
            ha="right",
            fontsize=7.2,
        )
        ax.set_yticks(
            np.arange(matrix.shape[0]), [str(int(x)) for x in matrix.index], fontsize=8
        )
        ax.set_xlabel("Cluster-associated gene")
        ax.set_ylabel("Transcriptional spot cluster")
        specimen_label(ax, sample, int(summary.loc[sample, "spots"]))
        bar = fig.colorbar(image, ax=ax, shrink=0.72, pad=0.018)
        bar.set_label("Mean log1p(CP10K)", fontsize=8, fontweight="bold")
        bar.ax.tick_params(labelsize=7)
    fig.subplots_adjust(
        left=0.055, right=0.985, top=0.975, bottom=0.085, wspace=0.28, hspace=0.42
    )
    return save(fig, "Figure04_Cluster_Associated_Marker_Expression")


def main():
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    outputs = []
    for maker in (figure02_umap, figure03_spatial_clusters, figure04_marker_expression):
        result = maker()
        outputs.extend(result)
        print(maker.__name__, *result, sep="\n  ")
    print(f"Completed 3 figures and {len(outputs)} exports.")


if __name__ == "__main__":
    main()
