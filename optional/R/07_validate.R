#!/usr/bin/env Rscript
##
## 07_validate.R -- independent checks on the estimated proportions.
##
##   Rscript 07_validate.R <params.json>
##
## Deconvolution fails quietly. Each check below can fail on its own, and a
## disagreement localises the problem rather than just flagging one.
##

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(BRETIGEA); library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 07_validate.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir <- params$output_dir
sfx     <- if (is.null(params$label) || !nzchar(params$label)) "" else paste0("_", params$label)
agg     <- params$aggregator
genes   <- params$check_genes
n_mark  <- params$n_marker
outf <- function(stem, ext) file.path(out_dir, sprintf("07_%s%s.%s", stem, sfx, ext))

expr  <- readRDS(file.path(out_dir, "02_expression.rds"))
props <- fread(file.path(out_dir, sprintf("05_cell_proportions%s.csv", sfx)))
props <- props[aggregator == agg]
if (!nrow(props)) stop("aggregator '", agg, "' not present in the proportions file")

## Match the matrix the dsa run actually used unless overridden.
dsa_sum  <- file.path(out_dir, sprintf("05_dsa_summary%s.json", sfx))
mat_name <- if (!is.null(params$expression)) {
  params$expression
} else if (file.exists(dsa_sum)) {
  fromJSON(dsa_sum)$expression
} else {
  "tpm"
}
mat <- expr[[mat_name]]
cat(sprintf("Validating the '%s' run against the %s matrix, %s aggregation.\n",
            if (nzchar(sfx)) sub("^_", "", sfx) else "unlabelled", mat_name, agg))

wide <- dcast(props, rids + con ~ cell_type, value.var = "proportion")
mat  <- mat[, wide$rids, drop = FALSE]

## ---- 1. direct marker expression -------------------------------------------
present <- intersect(genes, rownames(mat))
missing <- setdiff(genes, rownames(mat))
marker_check <- rbindlist(lapply(present, function(g) {
  x <- mat[g, ]
  data.table(gene = g,
             median_case = median(x[wide$con == sort(unique(wide$con))[1]]),
             median_ref  = median(x[wide$con == sort(unique(wide$con))[2]]))
}))
if (nrow(marker_check)) {
  marker_check[, ratio := median_case / median_ref]
  fwrite(marker_check, outf("marker_check", "csv"))
}

## ---- 2. BRETIGEA's SVD estimator -------------------------------------------
## Different mathematics, same markers. Output is surrogate variables on an
## arbitrary scale, so only the rank correlation across samples is meaningful.
cts <- setdiff(names(wide), c("rids", "con"))
svd_cmp <- tryCatch({
    sv <- NULL
  invisible(capture.output(
    sv <- BRETIGEA::brainCells(as.matrix(log2(mat + 1)), nMarker = n_mark,
                               species = "human", celltypes = cts)))
  rbindlist(lapply(intersect(cts, colnames(sv)), function(ct)
    data.table(cell_type = ct,
               spearman = suppressWarnings(
                 cor(wide[[ct]], sv[, ct], method = "spearman")))))
}, error = function(e) {
  cat("BRETIGEA::brainCells failed:", conditionMessage(e), "\n"); NULL
})
if (!is.null(svd_cmp)) fwrite(svd_cmp, outf("svd_comparison", "csv"))

## ---- 3. PC1 against the neuronal proportion --------------------------------
lg <- log2(mat + 1)
v  <- head(order(rowVars(lg), decreasing = TRUE), min(2000L, nrow(lg)))
pc <- prcomp(t(lg[v, ]), scale. = TRUE)
pct1 <- 100 * pc$sdev[1]^2 / sum(pc$sdev^2)
pc1_cmp <- rbindlist(lapply(cts, function(ct)
  data.table(cell_type = ct,
             spearman = suppressWarnings(
               cor(pc$x[, 1], wide[[ct]], method = "spearman")))))
pc1_cmp[, pc1_variance_pct := pct1]
fwrite(pc1_cmp, outf("pc1_comparison", "csv"))

## The three checks above are all internal: they compare this run against other
## quantities derived from this same cohort, so each can fail on its own and a
## disagreement localises the problem.
##
## There is deliberately no fourth check against literature proportions. A
## fold-change against a number from a different cohort, region and method is
## not a validation -- it reads as one while testing nothing. Comparisons to
## the literature belong in prose, where the provenance of each figure can be
## stated and checked.

## ---- plot ------------------------------------------------------------------
png(outf("validation", "png"), width = 1600, height = 800, res = 150)
par(mfrow = c(1, 2))
if ("neu" %in% cts) {
  plot(pc$x[, 1], wide$neu, pch = 19, col = factor(wide$con),
       xlab = sprintf("PC1 (%.1f%%)", pct1), ylab = "neuronal proportion",
       main = "PC1 vs neuronal proportion")
}
if (!is.null(svd_cmp) && nrow(svd_cmp)) {
  barplot(svd_cmp$spearman, names.arg = svd_cmp$cell_type, ylim = c(-1, 1),
          ylab = "Spearman, DSA vs BRETIGEA SVD", main = "Orthogonal estimator")
  abline(h = c(0, 0.8), lty = c(1, 2), col = c("black", "red"))
}
invisible(dev.off())

## ---- report ----------------------------------------------------------------
cat("\n== 1. Direct marker expression by group ==\n")
if (nrow(marker_check)) {
  print(marker_check[, .(gene, median_case = round(median_case, 1),
                         median_ref = round(median_ref, 1),
                         ratio = round(ratio, 3))], row.names = FALSE)
} else cat("none of the check genes are in the matrix\n")
if (length(missing)) cat(sprintf("not in matrix: %s\n", paste(missing, collapse = ", ")))
cat("If a marker is flat while its cell type shows a large estimated shift,\n")
cat("the estimate is not supported by the underlying expression.\n")

cat("\n== 2. BRETIGEA SVD estimator (Spearman across samples) ==\n")
if (!is.null(svd_cmp)) {
  print(svd_cmp[, .(cell_type, spearman = round(spearman, 3))], row.names = FALSE)
  cat("Above 0.8 is reassuring. Rank correlation only -- the scales differ.\n")
} else cat("unavailable\n")

cat(sprintf("\n== 3. PC1 (%.1f%% of variance) vs proportions ==\n", pct1))
print(pc1_cmp[, .(cell_type, spearman = round(spearman, 3))], row.names = FALSE)
cat("A strong neuronal correlation means two independent routes agree on the\n")
cat("dominant axis of the transcriptome.\n")

write_json(list(label = params$label, expression = mat_name, aggregator = agg,
                genes_missing = missing, marker_check = marker_check,
                svd = svd_cmp, pc1 = pc1_cmp, pc1_variance_pct = pct1),
           outf("validate_summary", "json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("\nWrote results to %s\n", out_dir))
