#!/usr/bin/env Rscript
##
## 03_sample_qc.R -- sex check, PCA and correlation outliers.
##
##   Rscript 03_sample_qc.R <params.json>
##
## Reports proposed exclusions; applies none. Dropping a sample changes the
## cohort, so that decision is the user's.
##

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 03_sample_qc.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir  <- params$output_dir
sexgenes <- params$sex_genes
n_var    <- params$n_variable
cor_sd   <- params$cor_sd

expr <- readRDS(file.path(out_dir, "02_expression.rds"))
tpm  <- expr$tpm
covs <- fread(file.path(out_dir, "03_covariates.csv"))
setkey(covs, rids)
covs <- covs[colnames(tpm)]

lg <- log2(tpm + 1)

## ---- sex check -------------------------------------------------------------
## No assumption about how msex is coded: samples are split by their own X/Y
## contrast, then cross-tabulated against the recorded value. A clean 2x2
## establishes the coding and mismatches fall off the diagonal.
x_genes <- intersect("XIST", rownames(lg))
y_genes <- intersect(setdiff(sexgenes, "XIST"), rownames(lg))
if (length(y_genes) == 0L) {
  warning("no Y-linked genes from --sex-genes are in the matrix; skipping sex check")
  sex_call <- rep(NA_character_, ncol(lg))
  x_score <- y_score <- rep(NA_real_, ncol(lg))
} else {
  y_score <- colMeans(lg[y_genes, , drop = FALSE])
  x_score <- if (length(x_genes)) lg[x_genes, ] else rep(0, ncol(lg))
  contrast <- y_score - x_score
  ## Two clear modes are expected; split at the midpoint of the gap.
  sorted <- sort(contrast)
  cut_at <- sorted[which.max(diff(sorted))] + max(diff(sorted)) / 2
  sex_call <- ifelse(contrast > cut_at, "Y_high", "X_high")
}

sex_tab <- if (all(is.na(sex_call))) NULL else table(sex_call, covs$msex)

## ---- PCA and correlation ---------------------------------------------------
## Named var_idx, not v: the driver loop below reuses its own loop variable, and
## a bare `v` there silently overwrote this one, so the summary reported the
## length of a column name instead of the number of genes.
var_idx <- head(order(rowVars(lg), decreasing = TRUE), min(n_var, nrow(lg)))
pc <- prcomp(t(lg[var_idx, ]), scale. = TRUE)
pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)

cm  <- cor(lg[var_idx, ], method = "spearman")
diag(cm) <- NA_real_
med_cor <- colMedians(cm, na.rm = TRUE)
cutoff  <- median(med_cor) - cor_sd * mad(med_cor)
outlier <- med_cor < cutoff

## ---- what the leading components track -------------------------------------
## The instrument for judging the annotation-based exclusion. PC1 correlating
## strongly with the removed fraction, or with library size, means the leading
## axis of the matrix is still library composition; correlating with condition
## instead is the outcome the exclusion was meant to expose.
drivers <- data.table(rids = colnames(tpm))
ef_file <- file.path(out_dir, "02_excluded_fraction.csv")
if (file.exists(ef_file)) {
  ef <- fread(ef_file); setkey(ef, rids)
  ef <- ef[colnames(tpm)]
  for (v in setdiff(names(ef), "rids")) drivers[[v]] <- as.numeric(ef[[v]])
}
drivers[, condition := as.integer(factor(covs$con))]

npc <- min(3L, ncol(pc$x))
pc_assoc <- rbindlist(lapply(setdiff(names(drivers), "rids"), function(v) {
  x <- drivers[[v]]
  if (sum(!is.na(x)) < 3L || length(unique(x[!is.na(x)])) < 2L) return(NULL)
  rbindlist(lapply(seq_len(npc), function(i) {
    ct <- suppressWarnings(cor.test(pc$x[, i], x, method = "spearman"))
    data.table(driver = v, pc = paste0("PC", i), pc_var_pct = pct[i],
               rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3))
  }))
}), fill = TRUE)
fwrite(pc_assoc, file.path(out_dir, "03_pc_associations.csv"))

qc <- data.table(
  rids       = colnames(tpm),
  con        = covs$con,
  msex       = covs$msex,
  x_score    = as.numeric(x_score),
  y_score    = as.numeric(y_score),
  sex_call   = sex_call,
  PC1        = pc$x[, 1],
  PC2        = pc$x[, 2],
  median_cor = med_cor,
  cor_outlier = outlier
)
fwrite(qc, file.path(out_dir, "03_sample_qc.csv"))

## ---- plots -----------------------------------------------------------------
png(file.path(out_dir, "03_pca.png"), width = 1600, height = 800, res = 150)
par(mfrow = c(1, 2))
grp <- factor(covs$con)
plot(pc$x[, 1], pc$x[, 2], col = grp, pch = 19,
     xlab = sprintf("PC1 (%.1f%%)", pct[1]), ylab = sprintf("PC2 (%.1f%%)", pct[2]),
     main = "PCA by condition")
legend("topright", levels(grp), col = seq_along(levels(grp)), pch = 19, bty = "n")
plot(sort(med_cor), pch = 19, ylab = "median Spearman to other samples",
     xlab = "sample (ranked)", main = "Sample correlation")
abline(h = cutoff, col = "red", lty = 2)
invisible(dev.off())

## ---- report ----------------------------------------------------------------
cat("\n== Sample QC ==\n")
cat(sprintf("samples            : %d\n", ncol(tpm)))
cat(sprintf("PC1/PC2 variance   : %.1f%% / %.1f%%\n", pct[1], pct[2]))

if (nrow(pc_assoc)) {
  cat("\n== What the leading components track (Spearman rho) ==\n")
  print(dcast(pc_assoc, driver ~ pc, value.var = "rho"), row.names = FALSE)
  cat("frac_excluded is the fraction step 2 removed by annotation. A high |rho|\n")
  cat("against PC1 means the exclusion did not go far enough: the leading axis\n")
  cat("still tracks how much was taken out rather than what remains.\n")
}

cat("\n== Sex check (expression call vs recorded msex) ==\n")
if (is.null(sex_tab)) {
  cat("skipped -- no Y-linked genes in the matrix\n")
} else {
  print(sex_tab)
  offdiag <- sum(sex_tab) - sum(apply(sex_tab, 1, max))
  cat(sprintf("\n%d sample(s) do not follow the dominant coding.\n", offdiag))
  if (offdiag > 0L) {
    mism <- qc[, .(rids, con, msex, sex_call, x_score = round(x_score, 2),
                   y_score = round(y_score, 2))]
    print(mism[order(sex_call, msex)], row.names = FALSE)
  }
}

cat(sprintf("\n== Correlation outliers (median cor < %.3f) ==\n", cutoff))
if (any(outlier)) {
  print(qc[cor_outlier == TRUE, .(rids, con, median_cor = round(median_cor, 3))],
        row.names = FALSE)
} else {
  cat("none\n")
}

write_json(list(
  n_samples      = ncol(tpm),
  n_variable     = length(var_idx),
  pc1_pct        = pct[1], pc2_pct = pct[2],
  cor_cutoff     = cutoff,
  n_cor_outliers = sum(outlier),
  sex_table      = if (is.null(sex_tab)) NULL else as.data.frame(sex_tab)
), file.path(out_dir, "03_sample_qc_summary.json"), auto_unbox = TRUE, pretty = TRUE)

cat("\nNo samples were excluded. Review 03_sample_qc.csv and decide.\n")
cat(sprintf("Wrote results to %s\n", out_dir))
