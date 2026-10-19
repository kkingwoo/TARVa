#!/usr/bin/env Rscript
## ============================================================================
## Figure: the magnitude of composition confounding.
##
## Five panels, all from this cohort's own results -- no reference values.
##   A  cell-type proportions by condition        (is there a confounder?)
##   B  compPC1 by condition                      (the axis everything uses)
##   C  per-site effect, before vs after          (did site findings move?)
##   D  per-gene effect, before vs after          (did gene findings move?)
##   E  significance, before vs after             (did anything cross alpha?)
##
## Panels C and D are separate because their coefficients are on different
## scales -- log-odds from the site models, rank units from the gene-level
## rank ANCOVA. Plotting them on shared axes would imply a comparability that
## does not exist. Panel E can combine them because p-values are scale-free.
##
## Usage: Rscript 11_confounding_plots.R <out_dir> [label]
## ============================================================================

suppressPackageStartupMessages({ library(data.table) })

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "."
sfx     <- if (length(args) >= 2 && nzchar(args[2])) paste0("_", args[2]) else ""

props <- fread(file.path(out_dir, sprintf("05_cell_proportions%s.csv", sfx)))
sites <- fread(file.path(out_dir, "09_principal_sites_adjusted.csv"))
genes <- fread(file.path(out_dir, "10_gene_level_adjusted.csv"))

ALPHA  <- 0.05
FLOOR  <- 1e-16          # logistf underflows to 0; -log10(0) is not plottable
COL    <- c("#0072B2", "#D55E00")   # reference, case -- colourblind-safe
PCH    <- c(presence = 21, proportion = 24)

use_agg <- props$aggregator[1]
p <- props[aggregator == use_agg]

grps <- sort(unique(p$con))
case <- grps[1]
ctrl <- setdiff(grps, case)
cols <- setNames(COL, c(ctrl, case))

## Display labels only. Every subset, model and colour lookup below keys on the
## database values in `grps`; this vector is consulted at draw time so the
## figure can say CN where the database says Control.
glab <- setNames(grps, grps)
glab[ctrl] <- "CN"

## compPC1, built exactly as steps 9 and 10 build it.
w <- dcast(p, rids + con ~ cell_type, value.var = "proportion")
P <- as.matrix(w[, setdiff(names(w), c("rids", "con")), with = FALSE])
rownames(P) <- w$rids
P[P <= 0] <- min(P[P > 0]) / 2
P   <- P / rowSums(P)
clr <- log(P) - rowMeans(log(P))
pca <- prcomp(clr, center = TRUE, scale. = FALSE)
w[, compPC1 := as.numeric(scale(pca$x[, 1]))]
ve1 <- pca$sdev[1]^2 / sum(pca$sdev^2) * 100
gate_p <- suppressWarnings(wilcox.test(compPC1 ~ con, data = w)$p.value)

## Points jittered onto every box: with 29 and 30 per group the raw spread is
## the evidence, and a box alone hides how far the groups overlap.
strip <- function(vals, grp, at, box_w = 0.55) {
  for (i in seq_along(at)) {
    v <- vals[grp == names(at)[i]]
    points(jitter(rep(at[i], length(v)), amount = box_w / 4), v,
           pch = 16, cex = 0.5, col = adjustcolor(cols[names(at)[i]], 0.55))
  }
}

## Labels for points that sit close together on a diagonal. Offsets every
## label to the right, then pushes them apart vertically until each clears the
## next, and recentres so the block does not drift off the cloud. A leader line
## keeps each label attached to its own point once it has moved.
label_points <- function(x, y, labs, cex = 0.62) {
  ord  <- order(y)
  x <- x[ord]; y <- y[ord]; labs <- labs[ord]
  gap  <- strheight("M", cex = cex) * 1.5
  ly   <- y
  for (i in seq_along(ly)[-1]) {
    if (ly[i] - ly[i - 1] < gap) ly[i] <- ly[i - 1] + gap
  }
  ly <- ly - (mean(ly) - mean(y))
  dx <- diff(par("usr")[1:2]) * 0.04
  segments(x, y, x + dx, ly, col = "grey60", lwd = 0.7)
  text(x + dx, ly, labels = labs, pos = 4, cex = cex, col = "#333333", xpd = NA)
}

open_dev <- function(path) {
  if (grepl("\\.png$", path)) png(path, width = 2400, height = 1500, res = 200)
  else pdf(path, width = 12, height = 7.5)
}

draw <- function() {
  layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 2, byrow = TRUE))
  par(mar = c(4.2, 4.4, 3.0, 1.0), cex.lab = 1.05)

  ## -- A. proportions by condition -------------------------------------------
  cts <- sort(unique(p$cell_type))
  at  <- seq_along(cts) * 3
  plot(NA, xlim = c(1, max(at) + 1.6), ylim = range(p$proportion) * c(0.9, 1.05),
       xaxt = "n", xlab = "", ylab = "estimated proportion",
       main = "A  Cell-type composition by condition")
  axis(1, at = at, labels = cts)
  for (i in seq_along(cts)) {
    for (j in seq_along(grps)) {
      v <- p[cell_type == cts[i] & con == grps[j], proportion]
      xj <- at[i] + (j - 1.5) * 0.75
      boxplot(v, at = xj, add = TRUE, boxwex = 0.6, axes = FALSE,
              outline = FALSE, border = cols[grps[j]], lwd = 1.6)
      points(jitter(rep(xj, length(v)), amount = 0.16), v, pch = 16, cex = 0.5,
             col = adjustcolor(cols[grps[j]], 0.5))
    }
  }
  legend("topright", legend = unname(glab[grps]), col = cols[grps], pch = 16, bty = "n")

  ## -- B. the adjustment axis ------------------------------------------------
  at2 <- setNames(c(1, 2), grps)
  plot(NA, xlim = c(0.4, 2.6), ylim = range(w$compPC1) * 1.15, xaxt = "n",
       xlab = "", ylab = "compPC1 (SD units)",
       main = sprintf("B  Composition axis\n%.0f%% of variance, p = %.3f",
                      ve1, gate_p))
  axis(1, at = at2, labels = unname(glab[names(at2)]))
  for (g in grps) {
    boxplot(w[con == g, compPC1], at = at2[g], add = TRUE, boxwex = 0.55,
            axes = FALSE, outline = FALSE, border = cols[g], lwd = 1.6)
  }
  strip(w$compPC1, w$con, at2)

  ## -- C. per-site effects ---------------------------------------------------
  s <- sites[!is.na(beta_unadj) & !is.na(beta_adj)]
  lim <- range(c(s$beta_unadj, s$beta_adj))
  plot(s$beta_unadj, s$beta_adj, xlim = lim, ylim = lim, asp = 1,
       xlab = "beta, unadjusted", ylab = "beta, composition-adjusted",
       main = sprintf("C  Per-site effects (n = %d)", nrow(s)),
       pch = PCH[s$finding], bg = adjustcolor("#666666", 0.6), col = "#333333")
  abline(0, 1, lty = 2, col = "grey40")
  legend("topleft", legend = names(PCH), pch = PCH, pt.bg = "#999999",
         bty = "n", cex = 0.9)

  ## -- D. per-gene effects ---------------------------------------------------
  ## Padded on the right so the labels have somewhere to go. asp = 1 is kept
  ## because the identity line only reads as "no change" at 45 degrees; it
  ## widens ylim to match rather than squashing the axes.
  lim <- range(c(genes$beta_unadj, genes$beta_adj))
  plot(genes$beta_unadj, genes$beta_adj,
       xlim = c(lim[1], lim[2] + diff(lim) * 0.34), ylim = lim, asp = 1,
       xlab = "beta, unadjusted (rank units)", ylab = "beta, adjusted",
       main = sprintf("D  Per-gene effects (n = %d)", nrow(genes)),
       pch = 21, bg = adjustcolor("#666666", 0.6), col = "#333333")
  abline(0, 1, lty = 2, col = "grey40")
  label_points(genes$beta_unadj, genes$beta_adj, genes$gene_name)

  ## -- E. significance -------------------------------------------------------
  ## Both analyses on one panel: p-values carry no units, so this is the only
  ## place the two are legitimately comparable.
  e <- rbind(
    sites[!is.na(p_unadj), .(pu = pmax(p_unadj, FLOOR), pa = pmax(p_adj, FLOOR),
                             src = "site")],
    genes[, .(pu = pmax(p_unadj, FLOOR), pa = pmax(p_adj, FLOOR), src = "gene")])
  e[, `:=`(x = -log10(pu), y = -log10(pa))]
  lim <- c(0, max(e$x, e$y) * 1.05)
  plot(e$x, e$y, xlim = lim, ylim = lim, asp = 1,
       xlab = expression(-log[10]~p~", unadjusted"),
       ylab = expression(-log[10]~p~", adjusted"),
       main = "E  Significance, before vs after",
       pch = ifelse(e$src == "site", 21, 22),
       bg = adjustcolor("#666666", 0.6), col = "#333333")
  abline(0, 1, lty = 2, col = "grey40")
  abline(h = -log10(ALPHA), v = -log10(ALPHA), lty = 3, col = "#D55E00")
  legend("bottomright", legend = c("site", "gene", expression(alpha == 0.05)),
         pch = c(21, 22, NA), lty = c(NA, NA, 3),
         col = c("#333333", "#333333", "#D55E00"), pt.bg = "#999999", bty = "n")
}

for (ext in c("pdf", "png")) {
  path <- file.path(out_dir, sprintf("11_confounding_figure.%s", ext))
  open_dev(path)
  draw()
  dev.off()
  cat(sprintf("Wrote %s\n", path))
}

## The numbers the figure is asserting, printed so the caption can quote them
## rather than being read off the panels.
s <- sites[!is.na(p_unadj)]
ssig <- s[p_unadj < ALPHA]
gsig <- genes[p_unadj < ALPHA]
cat(sprintf("\nComposition axis vs condition: p = %.3f (compPC1, %.1f%% of variance)\n",
            gate_p, ve1))
cat(sprintf("Sites: %d significant unadjusted, %d after; median |change| %.1f%%\n",
            nrow(ssig), sum(ssig$p_adj < ALPHA), median(abs(ssig$pct_change))))
cat(sprintf("Genes: %d significant unadjusted, %d after; median |change| %.1f%%\n",
            nrow(gsig), sum(gsig$p_adj < ALPHA), median(abs(gsig$pct_change))))
