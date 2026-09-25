# Week 5 Homework -- DESeq2 analysis script
#
# Course: Bioinformatics: From Multi-Omics Data to Discovery
# Week 5: Transcriptomics -- RNA-seq principles and DESeq2
#
# This script is written to run from an ordinary Rscript command, for example:
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" week5_deseq2_analysis.R
#
# It expects these input files in the same folder:
#   Week5_Homework_Count_Matrix.csv
#   Week5_Homework_Sample_Metadata.csv
#
# Required packages:
#   BiocManager::install(c("DESeq2", "apeglm"))
#   install.packages(c("tidyverse", "ggrepel", "pheatmap"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
})

# ----------------------------------------------------------------------
# 1. Work from the for_student folder
# ----------------------------------------------------------------------
work_dir <- "C:/Users/33846/Desktop/大三上/生物信息学/week5/Homework/for_student"
if (!dir.exists(work_dir)) {
  stop("Working folder not found: ", work_dir)
}
setwd(work_dir)

count_file    <- "Week5_Homework_Count_Matrix.csv"
metadata_file <- "Week5_Homework_Sample_Metadata.csv"

if (!file.exists(count_file))    stop("Count file missing: ", count_file)
if (!file.exists(metadata_file)) stop("Metadata file missing: ", metadata_file)

dir.create("outputs",  showWarnings = FALSE, recursive = TRUE)
dir.create("figures",  showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------
# 2. Import
# ----------------------------------------------------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

# ----------------------------------------------------------------------
# 3. Mandatory validation
# ----------------------------------------------------------------------
stopifnot(ncol(counts) == nrow(coldata))
stopifnot(identical(colnames(counts), rownames(coldata)))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))

coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch     <- factor(coldata$batch)

cat("=== Sample table (batch x condition) ===\n")
print(table(coldata$batch, coldata$condition))

cat("\n=== Library sizes (rounded) ===\n")
print(summary(round(colSums(counts))))

# ----------------------------------------------------------------------
# 4. Build DESeq2 object
# ----------------------------------------------------------------------
# Why include batch?
# The samples were processed in three batches (A, B, C), and each batch
# contains both control and treated samples. If we omitted batch, any
# systematic difference between batches could be absorbed into the
# condition effect, which would make the treated-vs-control comparison
# harder to interpret. The design ~ batch + condition compares treated
# and control samples within the same batch first, then combines the
# information across batches in a balanced way.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ batch + condition
)

# ----------------------------------------------------------------------
# 5. Pre-filter low-count genes
# ----------------------------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3
cat("\nGenes before filtering:", nrow(dds), "\n")
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# ----------------------------------------------------------------------
# 6. Fit model
# ----------------------------------------------------------------------
dds <- DESeq(dds)

coef_names <- resultsNames(dds)
cat("\n=== Coefficient names from resultsNames(dds) ===\n")
print(coef_names)

# Pick the treated-vs-control coefficient automatically.
# With design = ~ batch + condition and control as reference,
# DESeq2 names the coefficient condition_treated_vs_control.
target_coef <- coef_names[grepl("condition_treated_vs_control", coef_names)]

if (length(target_coef) == 0) {
  stop(
    "Could not find the treated-vs-control coefficient.\n",
    "resultsNames(dds) =", paste(coef_names, collapse = ", "), "\n",
    "Please inspect the output above and set target_coef manually."
  )
}

cat("\nUsing target coefficient:", target_coef, "\n")

# ----------------------------------------------------------------------
# 7. Extract results and apply shrinkage
# ----------------------------------------------------------------------
res <- results(
  dds,
  contrast = c("condition", "treated", "control"),
  alpha    = 0.05
)

res_shrunk <- lfcShrink(
  dds,
  coef = target_coef,
  type = "apeglm"
)

res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  mutate(
    baseMean         = as.numeric(baseMean),
    log2FoldChange   = as.numeric(log2FoldChange),
    lfcSE            = as.numeric(lfcSE),
    pvalue           = as.numeric(pvalue),
    padj             = as.numeric(padj),
    significant      = !is.na(padj) &
                       padj < 0.05 &
                       abs(log2FoldChange) >= 1,
    direction = case_when(
      significant & log2FoldChange > 0 ~ "Up in treated",
      significant & log2FoldChange < 0 ~ "Down in treated",
      TRUE                               ~ "Not significant"
    )
  ) |>
  arrange(padj)

write.csv(
  res_df,
  "outputs/week5_deseq2_results.csv",
  row.names = FALSE
)

sig_total <- sum(res_df$significant)
sig_up    <- sum(res_df$direction == "Up in treated")
sig_down  <- sum(res_df$direction == "Down in treated")

cat("\n=== Significant genes ===\n")
cat("Total (padj < 0.05 and |log2FC| >= 1):", sig_total, "\n")
cat("Up in treated:", sig_up, "\n")
cat("Down in treated:", sig_down, "\n")
cat("\nDirection table:\n")
print(table(res_df$direction))

# ----------------------------------------------------------------------
# 8. PCA (use varianceStabilizingTransformation directly, since the
#    filtered matrix is smaller than the default nsub for vst())
# ----------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vsd_mat <- assay(vsd)

pca_obj <- prcomp(t(vsd_mat), center = TRUE, scale. = FALSE)
pca_df <- as.data.frame(pca_obj$x)
pca_df$name     <- rownames(pca_df)
pca_df$condition <- coldata$condition[rownames(pca_df)]
pca_df$batch     <- coldata$batch[rownames(pca_df)]

percent_var <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = condition, shape = batch)
) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = Inf) +
  labs(
    title    = "Week 5 RNA-seq PCA",
    x        = paste0("PC1: ", percent_var[1], "% variance"),
    y        = paste0("PC2: ", percent_var[2], "% variance"),
    color    = "Condition",
    shape    = "Batch"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "figures/week5_pca.png",
  p_pca,
  width = 7, height = 5, dpi = 300
)

# ----------------------------------------------------------------------
# 9. Volcano plot
# ----------------------------------------------------------------------
plot_df <- res_df |>
  mutate(
    neg_log10_padj = -log10(pmax(padj, 1e-300))
  )

p_volcano <- ggplot(
  plot_df,
  aes(x = log2FoldChange, y = neg_log10_padj, color = direction)
) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Up in treated"      = "#C0392B",
      "Down in treated"    = "#2F6DB3",
      "Not significant"    = "grey70"
    )
  ) +
  labs(
    title    = "Differential expression: treated versus control",
    x        = "Shrunken log2 fold change",
    y        = "-log10 adjusted p value",
    color    = NULL
  ) +
  theme_bw(base_size = 12)

ggsave(
  "figures/week5_de_plot.png",
  p_volcano,
  width = 7, height = 5, dpi = 300
)

# ----------------------------------------------------------------------
# 10. Save reproducibility files
# ----------------------------------------------------------------------
saveRDS(dds, "outputs/week5_deseq2_object.rds")

capture.output(sessionInfo(), file = "outputs/session_info.txt")

# ----------------------------------------------------------------------
# 11. Copy finished deliverables into the week5 folder
# ----------------------------------------------------------------------
dest_dir <- "C:/Users/33846/Desktop/大三上/生物信息学/week5"
if (!dir.exists(dest_dir)) dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

copy_file <- function(src, dst_name) {
  dst <- file.path(dest_dir, dst_name)
  file.copy(src, dst, overwrite = TRUE)
  cat("Copied:", src, "->", dst, "\n")
}

copy_file("outputs/week5_deseq2_results.csv",      "week5_deseq2_results.csv")
copy_file("outputs/week5_deseq2_object.rds",       "week5_deseq2_object.rds")
copy_file("outputs/session_info.txt",              "session_info.txt")
copy_file("figures/week5_pca.png",                 "week5_pca.png")
copy_file("figures/week5_de_plot.png",             "week5_de_plot.png")

cat("\n=== Done. Deliverables are in:\n")
cat("  ", dest_dir, "\n")
