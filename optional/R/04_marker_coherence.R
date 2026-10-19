#!/usr/bin/env Rscript
##
## 04_marker_coherence.R -- do a cell type's markers co-vary across samples?
##
##   Rscript 04_marker_coherence.R <params.json>
##
## Computed on log2 (correlations only); DSA still consumes linear TPM.
##

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(BRETIGEA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 04_marker_coherence.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir <- params$output_dir
top_n   <- params$top_n
coh_min <- params$coherence_min

expr <- readRDS(file.path(out_dir, "02_expression.rds"))
lg   <- log2(expr$tpm + 1)

md  <- BRETIGEA::markers_df_human_brain
cts <- sort(unique(md$cell))

coherence <- rbindlist(lapply(cts, function(ct) {
  g <- intersect(head(md$markers[md$cell == ct], top_n), rownames(lg))
  if (length(g) < 3L) {
    return(data.table(cell_type = ct, n_markers = length(g),
                      median_r = NA_real_, q1_r = NA_real_, q3_r = NA_real_,
                      min_r = NA_real_, max_r = NA_real_))
  }
  cc <- cor(t(lg[g, ]), method = "pearson")
  vals <- cc[upper.tri(cc)]
  data.table(cell_type = ct,
             n_markers = length(g),
             median_r  = median(vals),
             q1_r      = quantile(vals, .25),
             q3_r      = quantile(vals, .75),
             min_r     = min(vals),
             max_r     = max(vals))
}))
coherence[, passes := !is.na(median_r) & median_r >= coh_min]

## ---- cross-type distinguishability -----------------------------------------
## Coherence is the diagonal: it says a type's markers co-vary. It cannot say
## they co-vary with THAT type rather than a neighbouring one -- OPC and
## oligodendrocyte marker sets grade into each other biologically, so an OPC set
## could be tight because it partly IS an oligodendrocyte set. A type is only
## separable if its diagonal clearly exceeds its off-diagonals; without that
## margin, nnls has no basis for splitting the shared signal.
pools <- lapply(cts, function(ct)
  intersect(head(md$markers[md$cell == ct], top_n), rownames(lg)))
names(pools) <- cts

cross <- matrix(NA_real_, length(cts), length(cts), dimnames = list(cts, cts))
for (i in seq_along(cts)) {
  for (j in seq_along(cts)) {
    if (j < i) { cross[i, j] <- cross[j, i]; next }
    ga <- pools[[i]]; gb <- pools[[j]]
    if (length(ga) < 3L || length(gb) < 3L) next
    cross[i, j] <- if (i == j) {
      cc <- cor(t(lg[ga, ])); median(cc[upper.tri(cc)])
    } else {
      median(cor(t(lg[ga, ]), t(lg[gb, ])))
    }
  }
}

coherence[, self_r    := diag(cross)[cell_type]]
coherence[, max_other := vapply(cell_type, function(ct)
                          max(cross[ct, setdiff(cts, ct)], na.rm = TRUE), numeric(1))]
coherence[, nearest   := vapply(cell_type, function(ct) {
                          o <- cross[ct, setdiff(cts, ct)]
                          setdiff(cts, ct)[which.max(o)] }, character(1))]
coherence[, margin    := self_r - max_other]

fwrite(coherence, file.path(out_dir, "04_marker_coherence.csv"))
fwrite(data.table(cell_type = rownames(cross), as.data.table(cross)),
       file.path(out_dir, "04_marker_cross_correlation.csv"))

png(file.path(out_dir, "04_coherence.png"), width = 1200, height = 800, res = 150)
bp <- barplot(coherence$median_r, names.arg = coherence$cell_type,
              xlab = "cell type",
              ylab = "median pairwise Pearson r", ylim = c(0, 1),
              col = ifelse(coherence$passes, "steelblue", "grey70"))
arrows(bp, coherence$q1_r, bp, coherence$q3_r, angle = 90, code = 3, length = 0.05)
## Median and IQR printed above each whisker so the figure carries its own
## values; this replaces the companion table in the results.
text(bp, coherence$q3_r + 0.07,
     labels = sprintf("%.3f\n(%.2f-%.2f)", coherence$median_r,
                      coherence$q1_r, coherence$q3_r),
     cex = 0.7)
abline(h = coh_min, col = "red", lty = 2)
invisible(dev.off())

cat("\n== Marker coherence ==\n")
print(coherence[, .(cell_type, n_markers,
                    median_r = round(median_r, 3),
                    IQR = paste0(round(q1_r, 2), "-", round(q3_r, 2)),
                    passes)], row.names = FALSE)

cat("\n== Cross-type marker correlation (diagonal = coherence) ==\n")
print(round(cross, 3))

cat("\n== Separation (self minus nearest other) ==\n")
print(coherence[, .(cell_type, self_r = round(self_r, 3),
                    nearest, max_other = round(max_other, 3),
                    margin = round(margin, 3))], row.names = FALSE)
tight <- coherence[!is.na(margin) & margin < 0.1, cell_type]
if (length(tight)) {
  cat(sprintf("\nCAUTION: %s have little margin over their nearest neighbour.\n",
              paste(tight, collapse = ", ")))
  cat("Coherent but possibly not separable -- nnls would split the shared signal\n")
  cat("arbitrarily. This is reported, not gated; use --cell-types on step 5 to act.\n")
}

passing <- coherence[passes == TRUE, cell_type]
failing <- coherence[passes == FALSE, cell_type]
cat(sprintf("\nPassing (>= %.2f): %s\n", coh_min, paste(passing, collapse = ", ")))
if (length(failing)) {
  cat(sprintf("Failing          : %s\n", paste(failing, collapse = ", ")))
  cat("Step 5 will exclude these unless --cell-types overrides it.\n")
}

write_json(list(coherence_min = coh_min,
                top_n         = top_n,
                passing       = passing,
                failing       = failing,
                low_margin    = tight,
                coherence     = coherence,
                cross         = as.data.frame(cross)),
           file.path(out_dir, "04_coherence_summary.json"),
           auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("\nWrote results to %s\n", out_dir))
