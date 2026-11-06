# Scientific visualisations for colorectal cancer metagenomics
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table", "ggplot2"))
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

figure_dir <- file.path(project_root(), "figures", "reproduced")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_research <- function() {
    theme_bw(base_size = 12) + theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black",
        face = "bold"), axis.title = element_text(colour = "black", face = "bold"), plot.title = element_text(colour = "black",
        face = "bold"), legend.title = element_text(colour = "black", face = "bold"), panel.grid.minor = element_blank())
}

save_figure <- function(plot, name, width = 10, height = 7) {
    ggsave(file.path(figure_dir, paste0(name, ".svg")), plot, width = width, height = height, units = "in",
        bg = "white")
    if (!capabilities("cairo") && !requireNamespace("ragg", quietly = TRUE)) {
        stop("A WebP-capable graphics device is required (for example, ragg)")
    }
    if (requireNamespace("ragg", quietly = TRUE)) {
        ggsave(file.path(figure_dir, paste0(name, ".webp")), plot, device = ragg::agg_webp, width = width,
            height = height, units = "in", dpi = 240, bg = "white")
    } else {
        ggsave(file.path(figure_dir, paste0(name, ".webp")), plot, device = "webp", width = width, height = height,
            units = "in", dpi = 240, bg = "white")
    }
}

species <- read_results("microbial_associations", "CRC_species_discovery_all.csv")
stopifnot(all(c("CLR_difference", "P_value", "FDR", "Scientific_name") %in% names(species)))
species[, Direction := fifelse(FDR < 0.05 & CLR_difference > 0, "Higher in CRC", fifelse(FDR < 0.05 &
    CLR_difference < 0, "Lower in CRC", "Not significant"))]
species[, logp := -log10(pmax(P_value, .Machine$double.xmin))]
volcano <- ggplot(species, aes(CLR_difference, logp, colour = Direction)) + geom_point(alpha = 0.75,
    size = 1.8) + geom_vline(xintercept = 0, linetype = 2) + scale_colour_manual(values = c(`Higher in CRC` = "#B64151",
    `Lower in CRC` = "#286C91", `Not significant` = "#A4A4A4")) + labs(title = "Microbial associations with colorectal cancer",
    x = "Adjusted CLR difference (CRC versus healthy)", y = expression(-log[10](italic(p))), colour = "Association") +
    theme_research()
save_figure(volcano, "microbial_association_volcano")

meta <- read_results("microbial_associations", "CRC_candidate_meta_analysis.csv")
stopifnot(all(c("Scientific_name", "Meta_estimate", "CI_lower", "CI_upper") %in% names(meta)))
meta[, Scientific_name := ifelse(is.na(Scientific_name), Species_ID, Scientific_name)]
meta[, Species := factor(Scientific_name, levels = Scientific_name[order(Meta_estimate)])]
forest <- ggplot(meta, aes(Meta_estimate, Species)) + geom_vline(xintercept = 0, linetype = 2) + geom_segment(aes(x = CI_lower,
    xend = CI_upper, yend = Species), linewidth = 0.65) + geom_point(size = 2.7, colour = "#286C91") +
    labs(title = "Cross-project microbial associations", x = "Random-effects CLR estimate (95% CI)",
        y = NULL) + theme_research()
save_figure(forest, "microbial_random_effects_forest", 11, 9)

cohorts <- read_results("microbial_associations", "CRC_candidate_cohort_effects.csv")
stopifnot(all(c("Species_ID", "Study", "Estimate") %in% names(cohorts)))
cohorts <- merge(cohorts, unique(meta[, .(Species_ID, Scientific_name)]), by = "Species_ID", all.x = TRUE)
cohorts[, Species := ifelse(is.na(Scientific_name), Species_ID, Scientific_name)]
cohort_heatmap <- ggplot(cohorts, aes(Study, Species, fill = Estimate)) + geom_tile(colour = "white",
    linewidth = 0.3) + scale_fill_gradient2(low = "#286C91", mid = "white", high = "#B64151", midpoint = 0) +
    labs(title = "Microbial association heterogeneity across projects", x = "Sequencing project", y = NULL,
        fill = "CLR effect") + theme_research() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(cohort_heatmap, "cross_project_microbial_heatmap", 12, 9)

functional <- read_results("functional_associations", "CRC_functional_meta_analysis.csv")
stopifnot(all(c("Module_ID", "Meta_estimate", "CI_lower", "CI_upper") %in% names(functional)))
if (!"Biological_implication" %in% names(functional)) functional[, Biological_implication := Module_ID]
functional[, Label := ifelse(is.na(Biological_implication) | !nzchar(Biological_implication), Module_ID,
    paste0(Module_ID, ": ", Biological_implication))]
functional <- functional[order(FDR, na.last = TRUE)][seq_len(min(.N, 20L))]
functional[, Label := factor(Label, levels = Label[order(Meta_estimate)])]
functional_plot <- ggplot(functional, aes(Meta_estimate, Label)) + geom_vline(xintercept = 0, linetype = 2) +
    geom_segment(aes(x = CI_lower, xend = CI_upper, yend = Label), linewidth = 0.65) + geom_point(size = 2.7,
    colour = "#286C91") + labs(title = "CRC-associated predicted microbial functions", x = "Random-effects CLR estimate (95% CI)",
    y = NULL) + theme_research()
save_figure(functional_plot, "functional_random_effects_forest", 13, 9)

clinical <- read_results("clinical_comparisons", "three_group_species_associations.csv")
stopifnot(all(c("Scientific_name", "Comparison", "CLR_difference", "FDR") %in% names(clinical)))
clinical[, Annotation := ifelse(!is.na(FDR) & FDR < 0.05, "*", "")]
clinical_plot <- ggplot(clinical, aes(Comparison, Scientific_name, fill = CLR_difference)) + geom_tile(colour = "white",
    linewidth = 0.3) + geom_text(aes(label = Annotation), colour = "black", fontface = "bold", size = 5) +
    scale_fill_gradient2(low = "#286C91", mid = "white", high = "#B64151", midpoint = 0) + labs(title = "Microbial associations across clinical groups",
    x = NULL, y = NULL, fill = "Adjusted CLR effect", caption = "* Benjamini-Hochberg FDR < 0.05") +
    theme_research() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_figure(clinical_plot, "healthy_adenoma_crc_heatmap", 11, 9)
