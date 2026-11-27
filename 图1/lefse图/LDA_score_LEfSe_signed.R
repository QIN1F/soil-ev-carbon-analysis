# ================================================================
# Signed LDA-score plot for species-level LEfSe-style analysis
# EV-CK: negative bars; EV-P: positive bars
# ================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ---- 1. Input and output paths ---------------------------------
otu_file   <- "upload/OTU-EV(2).csv"
group_file <- "upload/groupE(2).txt"
out_dir    <- "lefse_results/signed_direction_R"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Required package: ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2 first: install.packages('ggplot2')")
}
library(ggplot2)

# ---- 2. Read and validate the data ------------------------------
otu <- read.csv(otu_file, check.names = FALSE, stringsAsFactors = FALSE)
group <- read.delim(group_file, header = TRUE, sep = "\t",
                    check.names = FALSE, stringsAsFactors = FALSE)

stopifnot("Taxonomy" %in% colnames(otu))
stopifnot(all(c("sample", "group1") %in% colnames(group)))

sample_names <- as.character(group$sample)
otu_samples <- setdiff(colnames(otu), "Taxonomy")
if (!setequal(sample_names, otu_samples)) {
  stop("Sample names in the OTU table and group file do not match.")
}

# Reorder OTU columns to exactly match the group file.
otu <- otu[, c("Taxonomy", sample_names), drop = FALSE]
abundance <- as.matrix(otu[, sample_names, drop = FALSE])
mode(abundance) <- "numeric"
if (anyNA(abundance)) stop("The OTU table contains missing or non-numeric values.")
if (any(abundance < 0)) stop("Abundance values cannot be negative.")

taxa <- as.character(otu$Taxonomy)
group$group1 <- factor(as.character(group$group1),
                       levels = unique(as.character(group$group1)))
if (nlevels(group$group1) != 2) stop("Exactly two groups are required.")
if (any(table(group$group1) < 2)) stop("Each group needs at least two replicates.")

group_a <- levels(group$group1)[1]  # EV-CK in the supplied metadata
group_b <- levels(group$group1)[2]  # EV-P in the supplied metadata
idx_a <- which(group$group1 == group_a)
idx_b <- which(group$group1 == group_b)

# Remove species that are zero in all samples.
keep <- rowSums(abundance) > 0
abundance <- abundance[keep, , drop = FALSE]
taxa <- taxa[keep]

# Library-size normalization followed by log10(CPM + 1), matching the
# previous analysis and the supplied rechecked LEfSe result table.
sample_totals <- colSums(abundance)
if (any(sample_totals <= 0)) stop("At least one sample has zero total abundance.")
cpm <- sweep(abundance, 2, sample_totals, "/") * 1e6
log_cpm <- log10(cpm + 1)

# ---- 3. Statistical tests and LEfSe-style effect size ------------
safe_kruskal <- function(x, y) {
  tryCatch(kruskal.test(list(x, y))$p.value, error = function(e) 1)
}

safe_wilcox <- function(x, y) {
  tryCatch(wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value,
           error = function(e) 1)
}

# For two groups, calculate the signed one-feature LDA coefficient as:
# (mean_EV-P - mean_EV-CK) / pooled within-group variance.
lda_effect <- function(x, y) {
  delta <- mean(y) - mean(x)
  pooled_var <- (((length(x) - 1) * var(x) +
                   (length(y) - 1) * var(y)) /
                 (length(x) + length(y) - 2))
  denominator <- max(pooled_var, 1e-12)
  signed_coef <- delta / denominator
  score_log10 <- log10(max(abs(signed_coef), 1e-12))
  c(coef_signed = signed_coef, score_log10 = score_log10)
}

result_list <- lapply(seq_len(nrow(log_cpm)), function(i) {
  x <- log_cpm[i, idx_a]
  y <- log_cpm[i, idx_b]
  eff <- lda_effect(x, y)
  mean_a <- mean(cpm[i, idx_a])
  mean_b <- mean(cpm[i, idx_b])
  enriched <- if (mean_b > mean_a) group_b else group_a
  data.frame(
    Taxonomy = taxa[i],
    Taxon_label = gsub("_", " ", taxa[i], fixed = TRUE),
    KW_P = safe_kruskal(x, y),
    Wilcoxon_P = safe_wilcox(x, y),
    mean_abundance_group_a_CPM = mean_a,
    mean_abundance_group_b_CPM = mean_b,
    LDA_coefficient_signed = unname(eff["coef_signed"]),
    LDA_score_log10 = unname(eff["score_log10"]),
    Enriched_group = enriched,
    stringsAsFactors = FALSE
  )
})

result <- do.call(rbind, result_list)
result$KW_FDR_BH <- p.adjust(result$KW_P, method = "BH")
result$Wilcoxon_FDR_BH <- p.adjust(result$Wilcoxon_P, method = "BH")
result$LEfSe_raw_significant <- result$KW_P < 0.05 & result$Wilcoxon_P < 0.05

significant <- result[result$LEfSe_raw_significant, , drop = FALSE]
significant$LDA_score_signed <- ifelse(
  significant$Enriched_group == group_b,
  significant$LDA_score_log10,
  -significant$LDA_score_log10
)

# Select the top 30 by absolute LEfSe-style LDA score.
significant <- significant[order(-significant$LDA_score_log10,
                                 significant$KW_P,
                                 significant$Wilcoxon_P), , drop = FALSE]
top30 <- head(significant, 30)
top30$Rank <- seq_len(nrow(top30))
top30 <- top30[, c("Rank", setdiff(colnames(top30), "Rank"))]

write.csv(result, file.path(out_dir, "lefse_all_species_results_R.csv"), row.names = FALSE)
write.csv(top30, file.path(out_dir, "lefse_top30_results_R.csv"), row.names = FALSE)

# ---- 4. Signed LDA score plot -----------------------------------
# Sort by signed score so EV-CK-enriched species are below zero and EV-P-
# enriched species are above zero. Labels retain the s__ species prefix.
plot_df <- top30[order(top30$LDA_score_signed), , drop = FALSE]
plot_df$Taxon_label <- factor(plot_df$Taxon_label, levels = plot_df$Taxon_label)
plot_df$Enriched_group <- factor(plot_df$Enriched_group,
                                 levels = c(group_a, group_b))

pal <- c("EV-CK" = "#93C4C6", "EV-P" = "#E79C9D")

p <- ggplot(plot_df, aes(x = Taxon_label, y = LDA_score_signed,
                          fill = Enriched_group)) +
  geom_col(width = 0.78, colour = "#4A4A4A", linewidth = 0.25) +
  geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.35) +
  scale_fill_manual(values = pal, drop = FALSE,
                    labels = c("EV-CK enriched", "EV-P enriched")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(
    title = "LEfSe discriminant species (top 30)",
    subtitle = paste0(
      "Species-level analysis; positive values indicate ", group_b,
      " enrichment and negative values indicate ", group_a,
      " enrichment; n = 6 per group"
    ),
    x = "Species",
    y = "LDA score",
    fill = "Enriched group"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 9, colour = "#555555", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1,
                               colour = "#222222", size = 8.5),
    axis.text.y = element_text(colour = "#222222"),
    legend.position = c(0.88, 0.90),
    legend.background = element_blank(),
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 9),
    panel.grid.major.y = element_line(colour = "#D9D9D9", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 16, 12)
  )

ggsave(file.path(out_dir, "LDA_score_top30_signed_R.pdf"), p,
       width = 20, height = 11, units = "cm", device = cairo_pdf)
ggsave(file.path(out_dir, "LDA_score_top30_signed_R.png"), p,
       width = 20, height = 11, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(out_dir, "LDA_score_top30_signed_R.svg"), p,
       width = 20, height = 11, units = "cm", bg = "white")

cat("Analysis completed.\n")
cat("Species tested:", nrow(result), "\n")
cat("Raw KW + Wilcoxon significant:", nrow(significant), "\n")
cat("Top 30 exported:", nrow(top30), "\n")
cat("EV-CK enriched in top 30:", sum(top30$Enriched_group == group_a), "\n")
cat("EV-P enriched in top 30:", sum(top30$Enriched_group == group_b), "\n")
