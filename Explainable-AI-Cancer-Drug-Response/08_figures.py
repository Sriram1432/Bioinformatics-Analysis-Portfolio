"""Regenerate the fifteen-figure cancer pharmacogenomics atlas.

Reads validated result tables and writes publication-oriented WebP and SVG figures with explicit scientific limitations.

Run from the repository root with ``python 08_figures.py``. Source files belong in
``Data`` beside this script and generated outputs are written to ``Results``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------
# Final scientific figure atlas
# ----------------------------------------------------------------------------

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D


def find_project_root() -> Path:
    """Locate Project5 when run as a file or from an interactive notebook/cell."""
    candidates: list[Path] = []
    if "__file__" in globals():
        candidates.append(PROJECT_ROOT)
    cwd = Path.cwd().resolve()
    candidates.extend([cwd, *cwd.parents])

    for candidate in candidates:
        if (candidate / "Results").is_dir():
            return candidate
    raise FileNotFoundError(
        "Could not locate the Project5 root. Start Python inside the Project5 folder "
        "or update the fallback path in find_project_root()."
    )


ROOT = find_project_root()
RESULTS = ROOT / "Results"
OUT = ROOT / "Final" / "Figures"
OUT.mkdir(parents=True, exist_ok=True)

INK = "#111827"
MUTED = "#5B6678"
GRID = "#DCE3EA"
NAVY = "#17324D"
BLUE = "#3973B9"
TEAL = "#008C86"
CYAN = "#56B4C2"
CORAL = "#E56B55"
GOLD = "#D6A441"
PURPLE = "#7B61A8"
LIGHT = "#F4F7FA"
MODALITY = {
    "RNA": BLUE,
    "PRO": TEAL,
    "CNV": GOLD,
    "MUT": CORAL,
    "AllFour": TEAL,
    "Genomic": PURPLE,
    "CancerMean": GOLD,
    "TrainMean": "#A5ADBA",
}

sns.set_theme(style="white", context="talk")
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 18,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "axes.labelcolor": INK,
        "axes.edgecolor": "#AEB9C5",
        "axes.linewidth": 0.8,
        "xtick.color": INK,
        "ytick.color": INK,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
    }
)

captions: list[dict[str, str]] = []


def read(rel: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(RESULTS / rel, **kwargs)


def clean(ax: plt.Axes, grid: str | None = None) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.7, alpha=0.75)
    else:
        ax.grid(False)


def header(ax: plt.Axes, title: str, subtitle: str) -> None:
    # Keep exported figures visually clean; methodological detail belongs in the caption file.
    ax.set_title(title, loc="left", color=INK, pad=18)


def save(fig: plt.Figure, stem: str, title: str, caption: str, sources: str) -> None:
    fig.savefig(OUT / f"{stem}.webp", dpi=300, bbox_inches="tight", pad_inches=0.25)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    captions.append(
        {"figure": stem, "title": title, "caption": caption, "sources": sources}
    )


def f01_coverage() -> None:
    d = read("Phase1/gdsc_multiomics_coverage.csv")
    cols = ["RNA", "PROTEOMICS", "CNV", "MUTATION_RECORD"]
    for c in cols:
        d[c] = d[c].astype(str).str.lower().eq("true")
    order = ["RNA", "PROTEOMICS", "CNV", "MUTATION_RECORD"]
    patterns = (
        d.groupby(order)
        .size()
        .reset_index(name="models")
        .sort_values("models", ascending=False)
        .head(10)
    )
    availability = pd.DataFrame(
        {"modality": order, "models": [int(d[c].sum()) for c in order]}
    )
    labels = {
        "RNA": "RNA",
        "PROTEOMICS": "Proteome",
        "CNV": "CNV",
        "MUTATION_RECORD": "Mutation data",
    }

    fig = plt.figure(figsize=(11.8, 6.8))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.70, 1.70],
        height_ratios=[1.05, 0.78],
        wspace=0.52,
        hspace=0.08,
    )
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    y = np.arange(len(availability))
    ax1.barh(y, availability.models, color=[BLUE, TEAL, GOLD, CORAL], height=0.62)
    ax1.set_yticks(y, [labels[x] for x in availability.modality])
    ax1.invert_yaxis()
    ax1.set_xlim(0, len(d) * 1.08)
    ax1.set_xlabel("GDSC2 cell lines with data")
    for yi, v in zip(y, availability.models):
        ax1.text(
            v + 10, yi, f"{v:,}", va="center", fontweight="bold", fontsize=10, color=INK
        )
    clean(ax1, "x")

    patterns = patterns.reset_index(drop=True)
    px = np.arange(len(patterns))
    ax2.bar(px, patterns.models, color=NAVY, width=0.68)
    for xi, v in zip(px, patterns.models):
        ax2.text(
            xi,
            v + max(patterns.models) * 0.025,
            f"{v:,}",
            ha="center",
            va="bottom",
            fontsize=7.5,
            fontweight="bold",
        )
    ax2.set_ylabel("Cell lines")
    ax2.set_xticks([])
    clean(ax2, "y")
    for i, row in patterns.iterrows():
        active = np.flatnonzero([bool(row[c]) for c in order])
        if len(active) > 1:
            ax3.plot(
                [i, i], [active.min(), active.max()], color="#AAB4C0", lw=2, zorder=1
            )
        for j, c in enumerate(order):
            ax3.scatter(
                i,
                j,
                s=72,
                color=MODALITY.get(c.split("_")[0], TEAL)
                if bool(row[c])
                else "#D6DCE3",
                edgecolor="white",
                linewidth=0.6,
                zorder=2,
            )
    ax3.set_xlim(-0.6, len(patterns) - 0.4)
    ax3.set_ylim(len(order) - 0.5, -0.5)
    ax3.set_xticks(px, [str(i + 1) for i in px])
    ax3.set_xlabel("Availability combination")
    ax3.set_yticks(range(4), [labels[c] for c in order])
    ax3.tick_params(axis="y", labelsize=8.5)
    ax3.grid(False)
    ax3.spines[:].set_visible(False)
    header(
        ax1,
        "Multi-omics coverage across GDSC2 cell lines",
        f"{len(d):,} unique GDSC2 models; availability is assay presence, not complete feature-level data",
    )
    caption = (
        "Observed assay coverage and the ten most frequent availability combinations. Mutation coverage means at least one "
        "record/coverage indicator and must not be interpreted as mutation-positive status. Source: Phase 1 integration outputs."
    )
    save(
        fig,
        "01_multiomics_coverage",
        "Multi-omics coverage across GDSC2 cell lines",
        caption,
        "Results/Phase1/gdsc_multiomics_coverage.csv",
    )


def f02_heterogeneity() -> None:
    d = read("Phase2C/drug_between_cancer_kruskal_fdr.csv")
    for c in ["epsilon2", "q_fdr", "n_tested", "n_cancers"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["epsilon2"]).nlargest(22, "epsilon2").sort_values("epsilon2")
    fig, ax = plt.subplots(figsize=(10.8, 7.8))
    y = np.arange(len(d))
    strength = np.clip(-np.log10(d.q_fdr), 0, 30)
    ax.hlines(y, 0, d.epsilon2, color=GRID, linewidth=3)
    sc = ax.scatter(
        d.epsilon2,
        y,
        s=55 + 2.5 * np.sqrt(d.n_tested),
        c=strength,
        cmap="viridis",
        vmin=0,
        vmax=max(10, strength.max()),
        edgecolor="white",
        linewidth=0.8,
    )
    ax.set_yticks(y, d.drug_name)
    ax.set_xlabel("Kruskal–Wallis epsilon-squared effect size")
    ax.set_xlim(0, max(0.47, d.epsilon2.max() + 0.035))
    clean(ax, "x")
    header(
        ax,
        "Cancer lineage explains substantial response heterogeneity",
        "Top 22 of 295 drug IDs by omnibus effect size; tests compare eligible cancer groups within each drug",
    )
    cb = fig.colorbar(sc, ax=ax, fraction=0.028, pad=0.02)
    cb.set_label("−log₁₀(BH q)", fontsize=9)
    caption = (
        "Top cross-lineage GDSC2 AUC heterogeneity effects. Circle size reflects the number of tested cell lines; colour "
        "reflects BH-adjusted significance across drug-level omnibus tests. Kruskal–Wallis significance can reflect distributional "
        "differences beyond medians and does not establish clinical selectivity."
    )
    save(
        fig,
        "02_cross_cancer_heterogeneity",
        "Cancer lineage explains substantial response heterogeneity",
        caption,
        "Results/Phase2C/drug_between_cancer_kruskal_fdr.csv",
    )


def f03_enrichment() -> None:
    d = read("Phase6A/focus_lower_auc_enrichment.csv", dtype={"drug_id": str})
    for c in ["rank_effect", "q_fdr", "n_cancer"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    shared = (
        d.groupby(["drug_id", "drug_name"]).focus.nunique().sort_values(ascending=False)
    )
    ids = list(shared[shared >= 2].index)
    extra = (
        d.sort_values("rank_effect", ascending=False)
        .drop_duplicates("focus")
        .set_index(["drug_id", "drug_name"])
        .index.tolist()
    )
    ids = list(dict.fromkeys(ids + extra))[:18]
    z = (
        d.set_index(["drug_id", "drug_name"])
        .loc[d.set_index(["drug_id", "drug_name"]).index.isin(ids)]
        .reset_index()
    )
    z["label"] = z.drug_name + " [" + z.drug_id + "]"
    order = z.groupby("label").rank_effect.max().sort_values().index.tolist()
    cancers = ["Colorectal", "Pancreatic", "Breast", "NSCLC", "SCLC"]
    fig, ax = plt.subplots(figsize=(11.5, 8.2))
    for _, r in z.iterrows():
        x = cancers.index(r.focus)
        y = order.index(r.label)
        ax.scatter(
            x,
            y,
            s=80 + 900 * r.rank_effect**2,
            c=-np.log10(max(r.q_fdr, 1e-30)),
            cmap="magma",
            vmin=1.3,
            vmax=10,
            edgecolor="white",
            linewidth=0.8,
        )
        ax.text(
            x,
            y,
            f"{r.rank_effect:.2f}",
            ha="center",
            va="center",
            fontsize=6.8,
            color="white" if -np.log10(max(r.q_fdr, 1e-30)) < 6.5 else INK,
            fontweight="bold",
        )
    ax.set_xticks(range(5), cancers)
    ax.set_yticks(range(len(order)), order)
    ax.set_xlim(-0.55, 4.55)
    ax.set_ylim(-0.7, len(order) - 0.3)
    clean(ax)
    ax.set_facecolor(LIGHT)
    header(
        ax,
        "Cancer-associated lower-AUC enrichment",
        "Selected shared and high-effect contrasts; area = rank effect, colour = BH significance, label = effect",
    )
    cb = fig.colorbar(
        ScalarMappable(norm=Normalize(1.3, 10), cmap="magma"),
        ax=ax,
        fraction=0.026,
        pad=0.02,
    )
    cb.set_label("−log₁₀(BH q)", fontsize=9)
    caption = (
        "Within-drug cancer-versus-rest Mann–Whitney rank effects for selected Phase 6A contrasts meeting q<0.05, "
        "rank effect≥0.30 and lower median AUC. Blank cells mean the contrast was not selected, not necessarily no effect. "
        "Different GDSC drug IDs are retained separately; lower AUC is an experimental viability response, not clinical efficacy."
    )
    save(
        fig,
        "03_drug_cancer_enrichment",
        "Cancer-associated lower-AUC enrichment",
        caption,
        "Results/Phase6A/focus_lower_auc_enrichment.csv",
    )


def f04_within_variability() -> None:
    d = read("Phase2C/within_cancer_response_variability.csv")
    for c in ["n", "iqr_auc", "median_auc"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d.n >= 20].nlargest(16, "iqr_auc").sort_values("iqr_auc")
    cancer_short = {
        "B-Cell Non-Hodgkin's Lymphoma": "B-cell NHL",
        "Acute Myeloid Leukemia": "AML",
        "Small Cell Lung Carcinoma": "SCLC",
        "Non-Small Cell Lung Carcinoma": "NSCLC",
        "B-Lymphoblastic Leukemia": "B-ALL",
    }
    d["cancer_label"] = d.cancer_type.replace(cancer_short)
    d["label"] = d.DRUG_NAME + " - " + d.cancer_label
    fig, ax = plt.subplots(figsize=(11.2, 8.4))
    y = np.arange(len(d))
    norm = Normalize(d.median_auc.min(), d.median_auc.max())
    cmap = plt.get_cmap("crest")
    colors = cmap(norm(d.median_auc.to_numpy()))
    ax.barh(y, d.iqr_auc, color=colors, height=0.64, edgecolor="white", linewidth=0.8)
    for yi, r in enumerate(d.itertuples()):
        ax.text(
            r.iqr_auc + 0.009,
            yi,
            f"n={int(r.n)}",
            va="center",
            fontsize=7.5,
            color=MUTED,
        )
    sc = ScalarMappable(norm=norm, cmap=cmap)
    ax.set_yticks(y, d.label)
    ax.set_xlabel("Within-cancer AUC interquartile range")
    clean(ax, "x")
    header(
        ax,
        "Large within-lineage response variation remains",
        "Top eligible drug-cancer groups with n≥20; circle size = cell-line count, colour = median AUC",
    )
    cb = fig.colorbar(sc, ax=ax, fraction=0.028, pad=0.02)
    cb.set_label("Median fitted AUC", fontsize=9)
    caption = (
        "Largest observed within-cancer GDSC2 AUC interquartile ranges among groups containing at least 20 cell lines. "
        "These distributions motivate within-lineage biomarker modelling; they do not define clinical sensitive/resistant categories."
    )
    save(
        fig,
        "04_within_cancer_variability",
        "Large within-lineage response variation remains",
        caption,
        "Results/Phase2C/within_cancer_response_variability.csv",
    )


def f05_phase2_prediction() -> None:
    d = read("Phase2/selected_model_vs_baseline.csv")
    for c in ["MAE", "baseline_MAE", "MAE_improvement", "R2", "Spearman"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.sort_values("MAE_improvement")
    fig, ax = plt.subplots(figsize=(10.9, 7.2))
    y = np.arange(len(d))
    ax.hlines(y, d.MAE, d.baseline_MAE, color=GRID, lw=4)
    ax.scatter(
        d.baseline_MAE,
        y,
        s=100,
        color=GOLD,
        label="Cancer-mean baseline",
        edgecolor="white",
    )
    ax.scatter(
        d.MAE,
        y,
        s=100,
        color=TEAL,
        label="Validation-selected model",
        edgecolor="white",
    )
    ax.set_yticks(y, d.drug_name)
    ax.set_xlabel("Held-out test MAE (AUC units; lower is better)")
    clean(ax, "x")
    ax.legend(frameon=False, loc="lower right")
    header(
        ax,
        "RNA models usually outperform cancer-mean prediction",
        "Original Phase 2 test results; model family selected on validation data for each drug ID",
    )
    caption = (
        "Paired held-out test MAE for the original RNA-based model and a training-derived cancer-mean baseline across 15 "
        "high-coverage drug IDs. The test set was subsequently inspected during project development, so these estimates are historical "
        "exploratory evidence rather than a still-pristine confirmatory holdout."
    )
    save(
        fig,
        "05_rna_prediction_vs_baseline",
        "RNA models usually outperform cancer-mean prediction",
        caption,
        "Results/Phase2/selected_model_vs_baseline.csv",
    )


def f06_multiomics_benchmark() -> None:
    d = read("Phase5B/foldwise_modality_summary.csv")
    for c in ["mean_MAE", "median_R2"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    order = ["AllFour", "Genomic", "RNA", "CancerMean", "TrainMean"]
    d = d.set_index("modality").loc[order].reset_index()
    labels = {
        "AllFour": "All four omics",
        "Genomic": "RNA + mutation + CNV",
        "RNA": "RNA only",
        "CancerMean": "Cancer mean",
        "TrainMean": "Training mean",
    }
    fig, ax = plt.subplots(figsize=(10.6, 5.8))
    x = np.arange(len(d))
    colors = [TEAL, PURPLE, BLUE, GOLD, "#9AA4B2"]
    ax.scatter(
        x,
        d.mean_MAE,
        s=[360, 260, 260, 260, 260],
        color=colors,
        edgecolor="white",
        lw=1.5,
        zorder=3,
    )
    ax.plot(x, d.mean_MAE, color=GRID, lw=2, zorder=1)
    for xi, v in zip(x, d.mean_MAE):
        ax.text(xi, v + 0.00075, f"{v:.5f}", ha="center", fontweight="bold", fontsize=9)
    short_labels = {
        "AllFour": "All four\nomics",
        "Genomic": "RNA + mutation\n+ CNV",
        "RNA": "RNA\nonly",
        "CancerMean": "Cancer\nmean",
        "TrainMean": "Training\nmean",
    }
    ax.set_xticks(x, [short_labels[m] for m in d.modality], rotation=0, ha="center")
    ax.tick_params(axis="x", labelsize=8.5, pad=7)
    ax.set_ylim(0.0635, 0.0815)
    ax.set_ylabel("Mean fold MAE (AUC units; lower is better)")
    clean(ax, "y")
    header(
        ax,
        "Leakage-conscious multi-omics benchmarking",
        "569 original TRAIN lines; 15 preselected drugs; 3×4 patient-grouped CV; fold-wise selection and imputation",
    )
    caption = (
        "Mean cross-validation MAE from Phase 5B. Modalities use matched held-out samples and fixed XGBoost settings. "
        "All-four produces the lowest mean error, but the average gain over RNA is modest and the drugs were selected earlier in the project."
    )
    fig.subplots_adjust(bottom=0.16)
    save(
        fig,
        "06_multiomics_benchmark",
        "Leakage-conscious multi-omics benchmarking",
        caption,
        "Results/Phase5B/foldwise_modality_summary.csv",
    )


def f07_drug_gains() -> None:
    d = read("Phase5B/foldwise_drug_summary.csv")
    for c in ["gain_vs_rna", "better_folds"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.sort_values("gain_vs_rna")
    y = np.arange(len(d))
    colors = np.where(d.gain_vs_rna >= 0, TEAL, CORAL)
    fig, ax = plt.subplots(figsize=(10.8, 7.2))
    ax.axvline(0, color=INK, lw=1)
    ax.barh(
        y, d.gain_vs_rna, color=colors, height=0.62, edgecolor="white", linewidth=0.7
    )
    for yi, r in enumerate(d.itertuples()):
        ax.text(
            r.gain_vs_rna + 0.00022,
            yi,
            f"{r.gain_vs_rna:+.4f}  ({int(r.better_folds)}/12)",
            va="center",
            ha="left",
            fontsize=8.2,
        )
    ax.set_yticks(y, d.drug_name)
    ax.set_xlabel("MAE(RNA) − MAE(All four); positive favours multi-omics")
    clean(ax, "x")
    header(
        ax,
        "Multi-omics benefit varies materially by drug",
        "Mean paired fold gain; annotation gives folds favouring all-four out of 12",
    )
    caption = (
        "Drug-specific paired MAE differences from Phase 5B. All-four improves mean MAE for 13/15 drug IDs, with the largest "
        "gain for Nutlin-3a. Repeated folds overlap and are descriptive evaluations, not independent biological replicates."
    )
    save(
        fig,
        "07_drugwise_multiomics_gain",
        "Multi-omics benefit varies materially by drug",
        caption,
        "Results/Phase5B/foldwise_drug_summary.csv",
    )


def f08_nutlin_shap() -> None:
    d = read("Phase3B/shap_feature_importance.csv")
    d = d[(d.drug_name == "Nutlin-3a (-)") & (d.model == "AllFour")].copy()
    d["mean_abs_SHAP"] = pd.to_numeric(d.mean_abs_SHAP, errors="coerce")
    d = d.nlargest(18, "mean_abs_SHAP").sort_values("mean_abs_SHAP")
    fig, ax = plt.subplots(figsize=(10.6, 7))
    y = np.arange(len(d))
    colors = [MODALITY.get(x, "#8793A1") for x in d.block]
    ax.barh(
        y, d.mean_abs_SHAP, color=colors, height=0.66, edgecolor="white", linewidth=0.6
    )
    ax.set_yticks(y, d.feature)
    ax.set_xlabel("Mean absolute SHAP value (validation predictions)")
    clean(ax, "x")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=MODALITY[k],
            label=k,
            markersize=9,
        )
        for k in ["RNA", "PRO", "CNV", "MUT"]
    ]
    ax.legend(handles=handles, frameon=False, ncol=4, loc="lower right")
    header(
        ax,
        "Recorded TP53 events dominate Nutlin-3a attribution",
        "Top 18 all-four features in the reused exploratory validation cohort; SHAP magnitude is not causal effect",
    )
    caption = (
        "Mean absolute SHAP attribution for the Phase 3B Nutlin-3a all-four XGBoost model. MUT:TP53 is the dominant predictive "
        "feature, consistent with MDM2–p53 biology. Mutation events are recorded-summary indicators and do not establish functional "
        "loss-of-function status; validation reuse and limited resampling make this exploratory."
    )
    save(
        fig,
        "08_nutlin_tp53_shap",
        "Recorded TP53 events dominate Nutlin-3a attribution",
        caption,
        "Results/Phase3B/shap_feature_importance.csv",
    )


def f09_slfn11() -> None:
    d = read("Phase4C/paired_adjusted_association_evidence.csv")
    d = d[
        (d.gene == "SLFN11") & d.drug.isin(["Camptothecin", "Irinotecan", "Topotecan"])
    ].copy()
    # Multiple IDs can share a displayed drug name; choose the strongest cross-platform-supported pair per name.
    for c in [
        "adjusted_rho_GDSC2",
        "adjusted_rho_PRISM",
        "q_adjusted_GDSC2",
        "q_adjusted_PRISM",
    ]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["joint_q"] = d[["q_adjusted_GDSC2", "q_adjusted_PRISM"]].max(axis=1)
    d = (
        d.sort_values("joint_q")
        .drop_duplicates("drug")
        .sort_values("adjusted_rho_GDSC2")
    )
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(10.4, 5.5))
    for yi, r in enumerate(d.itertuples()):
        ax.plot(
            [r.adjusted_rho_GDSC2, r.adjusted_rho_PRISM], [yi, yi], color=GRID, lw=4
        )
    ax.scatter(
        d.adjusted_rho_GDSC2, y, s=180, color=BLUE, label="GDSC2", edgecolor="white"
    )
    ax.scatter(
        d.adjusted_rho_PRISM, y, s=180, color=CORAL, label="PRISM", edgecolor="white"
    )
    for yi, r in enumerate(d.itertuples()):
        ax.text(
            r.adjusted_rho_GDSC2 - 0.012,
            yi - 0.18,
            f"{r.adjusted_rho_GDSC2:.2f}",
            ha="center",
            color=BLUE,
            fontweight="bold",
            fontsize=9,
        )
        ax.text(
            r.adjusted_rho_PRISM + 0.012,
            yi - 0.18,
            f"{r.adjusted_rho_PRISM:.2f}",
            ha="center",
            color=CORAL,
            fontweight="bold",
            fontsize=9,
        )
    ax.axvline(0, color=INK, lw=0.8, ls="--")
    ax.set_yticks(y, d.drug)
    ax.set_xlim(-0.55, 0.04)
    ax.set_xlabel("Cancer- and RNA-source-adjusted partial rank association (ρ)")
    clean(ax, "x")
    ax.legend(frameon=False, loc="lower left", ncol=2)
    header(
        ax,
        "SLFN11 shows concordant TOP1-inhibitor associations",
        "Higher expression associates with lower fitted AUC on both platforms; all displayed pairs pass BH q<0.05",
    )
    caption = (
        "Cancer- and RNA-source-adjusted partial-rank associations between SLFN11 expression and drug AUC in overlapping GDSC2/PRISM "
        "cell lines. Compound matching remains name-provisional and the same models contribute to both assays; this is cross-platform "
        "concordance, not independent unseen-line validation or proof of mechanism."
    )
    save(
        fig,
        "09_slfn11_crossplatform",
        "SLFN11 shows concordant TOP1-inhibitor associations",
        caption,
        "Results/Phase4C/paired_adjusted_association_evidence.csv",
    )


def f10_prism_concordance() -> None:
    d = read(
        "Phase4A/cross_platform_shared_line_concordance.csv", dtype={"drug_id": str}
    )
    for c in ["spearman_auc", "q_fdr", "n_shared"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d.n_shared >= 50].nlargest(18, "spearman_auc").sort_values("spearman_auc")
    prism_suffix = d.broad_id.astype(str).str.split("-").str[-3:].str.join("-")
    d["label"] = d.gdsc_name + " [" + d.drug_id + "] · " + prism_suffix
    fig, ax = plt.subplots(figsize=(10.8, 7.8))
    y = np.arange(len(d))
    strength = np.clip(-np.log10(d.q_fdr), 0, 20)
    norm = Normalize(strength.min(), strength.max())
    cmap = plt.get_cmap("viridis")
    ax.barh(
        y,
        d.spearman_auc,
        color=cmap(norm(strength)),
        height=0.68,
        edgecolor="white",
        linewidth=0.7,
    )
    sc = ScalarMappable(norm=norm, cmap=cmap)
    ax.set_yticks(y, d.label)
    ax.set_xlabel("Spearman correlation of GDSC2 and PRISM AUC on shared lines")
    clean(ax, "x")
    header(
        ax,
        "Cross-platform response concordance is compound-dependent",
        "Top provisional name-matched pairs with ≥50 shared cell lines; circle size = shared n",
    )
    cb = fig.colorbar(sc, ax=ax, fraction=0.028, pad=0.02)
    cb.set_label("−log₁₀(BH q)", fontsize=9)
    caption = (
        "Highest shared-cell-line GDSC2–PRISM AUC correlations among provisional exact-name compound matches. Multiple GDSC drug IDs may "
        "share a name. None of these pairs is structure-certified in the current project, and overlap of cell lines prevents interpretation "
        "as independent external validation."
    )
    save(
        fig,
        "10_gdsc_prism_concordance",
        "Cross-platform response concordance is compound-dependent",
        caption,
        "Results/Phase4A/cross_platform_shared_line_concordance.csv",
    )


def f11_mapk_enrichment() -> None:
    d = read("Phase6A/focus_drug_cancer_enrichment_all.csv", dtype={"drug_id": str})
    ids = ["1014", "1060", "1062", "1372", "1564", "1736", "1908"]
    d = d[d.drug_id.isin(ids) & d.focus.isin(["Colorectal", "Pancreatic"])].copy()
    for c in ["rank_effect", "q_fdr", "n_cancer", "cancer_median_AUC"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["label"] = d.drug_name + " [" + d.drug_id + "]"
    order = d.drop_duplicates("drug_id").sort_values("drug_id").label.tolist()
    fig, ax = plt.subplots(figsize=(12.2, 7.8))
    for _, r in d.iterrows():
        x = ["Colorectal", "Pancreatic"].index(r.focus)
        y = order.index(r.label)
        sig = -np.log10(max(r.q_fdr, 1e-30))
        ax.scatter(
            x,
            y,
            s=70 + 900 * abs(r.rank_effect) ** 2,
            c=r.rank_effect,
            cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-0.2, vcenter=0, vmax=0.65),
            edgecolor=INK if r.q_fdr < 0.05 else "white",
            lw=1.5,
        )
        ax.text(
            x,
            y,
            f"{r.rank_effect:+.2f}\n(n={int(r.n_cancer)})",
            ha="center",
            va="center",
            fontsize=6.4,
            fontweight="bold",
            color="white" if abs(r.rank_effect) > 0.37 else INK,
        )
    ax.set_xticks([0, 1], ["Colorectal", "Pancreatic"])
    ax.set_yticks(range(len(order)), order)
    ax.set_xlim(-0.55, 1.55)
    clean(ax)
    ax.set_facecolor(LIGHT)
    header(
        ax,
        "MEK/ERK lower-AUC enrichment spans CRC and pancreatic lines",
        "Within-drug cancer-versus-rest rank effect; outlined circles pass global Phase 6A BH q<0.05",
    )
    caption = (
        "Colorectal and pancreatic rank-enrichment effects for seven selected MEK/ERK drug IDs. Positive values mean relatively lower AUC "
        "than other cancer lines for that drug. Pancreatic significance can coexist with high absolute AUC and does not establish a clinically "
        "actionable dependency; two rows are distinct selumetinib screening IDs."
    )
    save(
        fig,
        "11_crc_pancreatic_mapk_enrichment",
        "MEK/ERK lower-AUC enrichment spans CRC and pancreatic lines",
        caption,
        "Results/Phase6A/focus_drug_cancer_enrichment_all.csv",
    )


def f12_mapk_molecular() -> None:
    d = read(
        "Phase6C/all_within_cancer_molecular_associations.csv", dtype={"drug_id": str}
    )
    for c in ["effect", "q_fdr", "n"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d.cancer.isin(["Colorectal Carcinoma", "Pancreatic Carcinoma"])].copy()
    selected = d[(d.q_fdr < 0.10) & (d.effect.abs() >= 0.30)].copy()
    selected["feature_short"] = selected.feature.str.replace(
        r"^(RNA|PRO|CNV|MUT):", "", regex=True
    )
    selected["feature_short"] = selected.feature_short.replace(
        {"MAPK_feedback_score": "MAPK RNA score"}
    )
    selected["label"] = selected.drug + " [" + selected.drug_id + "]"
    features = (
        selected.groupby("feature_short")
        .effect.apply(lambda x: x.abs().max())
        .sort_values(ascending=False)
        .head(14)
        .index
    )
    selected = selected[selected.feature_short.isin(features)]
    drugs = (
        selected.groupby("label")
        .effect.apply(lambda x: x.abs().max())
        .sort_values()
        .index.tolist()
    )
    features = list(features[::-1])
    fig, ax = plt.subplots(figsize=(12, 7.2))
    for _, r in selected.iterrows():
        x = features.index(r.feature_short)
        y = drugs.index(r.label)
        sig = -np.log10(max(r.q_fdr, 1e-30))
        ax.scatter(
            x,
            y,
            s=60 + 650 * abs(r.effect) ** 2,
            c=r.effect,
            cmap="coolwarm",
            norm=TwoSlopeNorm(vmin=-0.7, vcenter=0, vmax=0.7),
            edgecolor=INK if r.q_fdr < 0.05 else "white",
            lw=1.3,
        )
    ax.set_xticks(range(len(features)), features, rotation=34, ha="right")
    ax.set_yticks(range(len(drugs)), drugs)
    clean(ax)
    fig.subplots_adjust(bottom=0.32)
    ax.set_facecolor(LIGHT)
    header(
        ax,
        "CRC MAPK associations are selective, not universal",
        "Target-cancer associations with |effect|≥0.30 and global BH q<0.10; outline denotes q<0.05",
    )
    cb = fig.colorbar(
        ScalarMappable(
            norm=TwoSlopeNorm(vmin=-0.7, vcenter=0, vmax=0.7), cmap="coolwarm"
        ),
        ax=ax,
        fraction=0.026,
        pad=0.02,
    )
    cb.set_label("Association effect", fontsize=9)
    caption = (
        "Selected within-cancer RNA/CNV/mutation associations for colorectal and pancreatic models. Continuous effects are Spearman ρ; "
        "mutation effects are median-AUC contrasts and are not directly comparable in scale. BH correction was applied across all 1,920 "
        "Phase 6C tests. No pancreatic association passed q<0.05."
    )
    save(
        fig,
        "12_mapk_molecular_associations",
        "CRC MAPK associations are selective, not universal",
        caption,
        "Results/Phase6C/all_within_cancer_molecular_associations.csv",
    )


def f13_mapk_signature_assoc() -> None:
    d = read(
        "Phase6C/all_within_cancer_molecular_associations.csv", dtype={"drug_id": str}
    )
    d = d[
        (d.feature == "RNA:MAPK_feedback_score")
        & d.cancer.isin(["Colorectal Carcinoma", "Pancreatic Carcinoma"])
    ].copy()
    for c in ["effect", "q_fdr", "n"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["label"] = d.drug + " [" + d.drug_id + "]"
    order = d[d.cancer == "Colorectal Carcinoma"].sort_values("effect").label.tolist()
    fig, ax = plt.subplots(figsize=(12.2, 7.9))
    for _, r in d.iterrows():
        x = 0 if r.cancer.startswith("Colorectal") else 1
        y = order.index(r.label)
        ax.scatter(
            x,
            y,
            s=100 + 900 * abs(r.effect) ** 2,
            c=r.effect,
            cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-0.65, vcenter=0, vmax=0.65),
            edgecolor=INK if r.q_fdr < 0.05 else "white",
            lw=1.6,
        )
        ax.text(
            x,
            y,
            f"{r.effect:+.2f}\nq={r.q_fdr:.2g}",
            ha="center",
            va="center",
            fontsize=6.2,
            fontweight="bold",
            color="white" if abs(r.effect) > 0.4 else INK,
        )
    ax.set_xticks([0, 1], ["Colorectal", "Pancreatic"])
    ax.set_yticks(range(len(order)), order)
    ax.set_xlim(-0.55, 1.55)
    clean(ax)
    ax.set_facecolor(LIGHT)
    header(
        ax,
        "MAPK feedback-score association is CRC-weighted",
        "Six-gene RNA score versus AUC; area = |Spearman ρ|; outline = global BH q<0.05",
    )
    caption = (
        "Within-cancer association of the source-standardised six-gene transcriptional feedback score with GDSC2 AUC. Three colorectal "
        "associations pass the global 1,920-test BH correction; none pass in pancreatic models. The score is an RNA feedback proxy and must "
        "not be labelled direct MAPK or phospho-ERK activity."
    )
    save(
        fig,
        "13_mapk_signature_associations",
        "MAPK feedback-score association is CRC-weighted",
        caption,
        "Results/Phase6C/all_within_cancer_molecular_associations.csv",
    )


def f14_signature_prediction() -> None:
    d = read("Phase6D/per_cancer_incremental_gains.csv", dtype={"drug_id": str})
    d["gain_signature_beyond_mutation"] = pd.to_numeric(
        d.gain_signature_beyond_mutation, errors="coerce"
    )
    wide = d.pivot(
        index=["drug_id", "drug"],
        columns="cancer_type",
        values="gain_signature_beyond_mutation",
    ).reset_index()
    wide = wide.sort_values("Colorectal Carcinoma").reset_index(drop=True)
    y = np.arange(len(wide))
    fig, ax = plt.subplots(figsize=(10.8, 6.6))
    ax.axvline(0, color=INK, lw=1)
    for yi, r in wide.iterrows():
        ax.plot(
            [r["Colorectal Carcinoma"], r["Pancreatic Carcinoma"]],
            [yi, yi],
            color=GRID,
            lw=4,
        )
    ax.scatter(
        wide["Colorectal Carcinoma"],
        y,
        s=140,
        color=TEAL,
        label="Colorectal",
        edgecolor="white",
    )
    ax.scatter(
        wide["Pancreatic Carcinoma"],
        y,
        s=140,
        color=CORAL,
        label="Pancreatic",
        edgecolor="white",
    )
    clean_name = wide.apply(
        lambda row: re.sub(rf"\s*\({re.escape(str(row.drug_id))}\)\s*$", "", str(row.drug)),
        axis=1,
    )
    ax.set_yticks(y, clean_name + " [" + wide.drug_id + "]")
    ax.set_xlabel("MAE(mutation model) − MAE(mutation + signature)")
    clean(ax, "x")
    ax.legend(frameon=False, loc="lower right")
    header(
        ax,
        "MAPK signature prediction depends on cancer context",
        "Positive gain favours adding the score; repeated grouped CV on 58 selected TRAIN lines",
    )
    caption = (
        "Cancer-specific out-of-fold MAE gain from adding the MAPK RNA score to cancer, RNA-source and recorded-mutation predictors. "
        "All seven drug IDs improve in colorectal lines, whereas six of seven worsen in pancreatic lines. Subgroups contain approximately "
        "34–36 and 19–22 lines; drugs/signature were selected after earlier exploration, so this is not independent confirmation."
    )
    save(
        fig,
        "14_mapk_signature_prediction_context",
        "MAPK signature prediction depends on cancer context",
        caption,
        "Results/Phase6D/per_cancer_incremental_gains.csv",
    )


def f15_ppi() -> None:
    d = read("Phase6B/ppi_vs_multiomics_fold_scores.csv")
    d.MAE = pd.to_numeric(d.MAE, errors="coerce")
    order = ["PD0325901", "Refametinib", "SCH772984"]
    models = ["AllFour", "PPI", "ShuffledPPI"]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    x = np.arange(len(models))
    drug_colors = {"PD0325901": BLUE, "Refametinib": CORAL, "SCH772984": TEAL}
    for drug in order:
        z = d[d.drug_name == drug].groupby("model").MAE.mean().reindex(models)
        c = drug_colors[drug]
        ax.plot(x, z.values, lw=2.4, alpha=0.86, color=c)
        ax.scatter(
            x, z.values, s=125, color=c, edgecolor="white", linewidth=1.1, zorder=3
        )
        ax.text(2.08, z.iloc[-1], drug, color=c, fontsize=9, weight="bold", va="center")
    ax.set_xticks(x, ["All four", "PPI informed", "Shuffled labels"])
    ax.set_xlim(-0.10, 2.62)
    ax.set_ylabel("Mean fold MAE (lower is better)")
    clean(ax, "y")
    header(
        ax,
        "STRING neighbourhood features add no consistent value",
        "2×4 patient-grouped CV; same fold-selected omics features and held-out lines",
    )
    caption = (
        "Mean Phase 6B fold MAE for standard all-four omics, added STRING v12.5 physical-neighbourhood summaries, and a shuffled-node-label "
        "control. Real PPI features are marginally better for one drug and worse for two. The shuffled control is not degree-preserving, so "
        "the result supports no consistent predictive benefit, not a general conclusion about network topology."
    )
    save(
        fig,
        "15_ppi_model_evaluation",
        "STRING neighbourhood features add no consistent value",
        caption,
        "Results/Phase6B/ppi_vs_multiomics_fold_scores.csv",
    )


def write_caption_index() -> None:
    frame = pd.DataFrame(captions)
    frame.to_csv(
        RESULTS / "Scientific_Audit" / "PROJECT5_FIGURE_CAPTIONS.csv", index=False
    )
    lines = ["# Project 5 figure captions", ""]
    for r in captions:
        lines += [
            f"## {r['figure']} — {r['title']}",
            "",
            r["caption"],
            "",
            f"**Data:** `{r['sources']}`",
            "",
        ]
    (RESULTS / "Scientific_Audit" / "PROJECT5_FIGURE_CAPTIONS.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    builders = [
        f01_coverage,
        f02_heterogeneity,
        f03_enrichment,
        f04_within_variability,
        f05_phase2_prediction,
        f06_multiomics_benchmark,
        f07_drug_gains,
        f08_nutlin_shap,
        f09_slfn11,
        f10_prism_concordance,
        f11_mapk_enrichment,
        f12_mapk_molecular,
        f13_mapk_signature_assoc,
        f14_signature_prediction,
        f15_ppi,
    ]
    for fn in builders:
        fn()
        print("built", fn.__name__)
    write_caption_index()
    print(f"Generated {len(captions)} standalone figures in {OUT}")


if __name__ == "__main__":
    main()
