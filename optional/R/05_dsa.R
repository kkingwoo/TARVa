#!/usr/bin/env Rscript
##
## 05_dsa.R -- Digital Sorting Algorithm proportions with Pottier's marker
## stability procedure.
##
##   Rscript 05_dsa.R <params.json>
##
## Reproduces CellMix's .DSAproportions(). The DSA package is GitHub-only and
## CellMix is off CRAN with a heavy dependency chain; the algorithm itself is
## six lines, so implementing it directly is the more reproducible route.
##

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(nnls); library(BRETIGEA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 05_dsa.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir  <- params$output_dir
top_n    <- params$top_n
n_mark   <- params$n_markers
n_iter   <- params$n_iter
seed     <- params$seed
aggs     <- params$aggregators
stab_min <- params$stability_min

## Optional suffix so alternative cell-type sets can coexist in one directory.
sfx <- if (is.null(params$label) || !nzchar(params$label)) "" else paste0("_", params$label)
outf <- function(stem, ext) file.path(out_dir, sprintf("05_%s%s.%s", stem, sfx, ext))

## ---- the algorithm ---------------------------------------------------------
## O is cell types x samples: the aggregated marker expression per type.
## Solve for a non-negative scaling g with  sum_k g_k * O_kj ~= 1  for every
## sample j, then p = g_k * O_kj is that type's share, renormalised to sum to 1.
## `O * g` scales row k by g[k] -- equivalent to diag(g) %*% O without the
## edge case where diag() of a length-1 vector builds an identity matrix.
dsa_proportions <- function(expr_lin, marker_list,
                            agg = c("median", "mean", "scaled_mean"),
                            renormalise = TRUE, prescaled = FALSE) {
  agg <- match.arg(agg)
  if (agg == "scaled_mean" && !prescaled) expr_lin <- expr_lin / rowMeans(expr_lin)
  f <- switch(agg,
              mean = colMeans, scaled_mean = colMeans,
              median = function(m) apply(m, 2, median))
  O <- t(vapply(marker_list,
                function(g) f(expr_lin[g, , drop = FALSE]),
                numeric(ncol(expr_lin))))
  rownames(O) <- names(marker_list); colnames(O) <- colnames(expr_lin)
  g <- coef(nnls(t(O), rep(1, ncol(O))))
  names(g) <- rownames(O)
  p <- t(O * g)
  dimnames(p) <- list(colnames(O), rownames(O))
  if (renormalise) p <- p / rowSums(p)
  ## nnls can set a coefficient to exactly zero when that cell type's marker
  ## profile is close to a linear combination of the others -- the algebraic
  ## form of poor separation (see step 4's margins). The type's proportion is
  ## then identically zero in every sample, which is degenerate rather than an
  ## estimate, so the caller needs to know it happened.
  attr(p, "zeroed") <- names(g)[g <= 0]
  p
}

## ---- inputs ----------------------------------------------------------------
expr     <- readRDS(file.path(out_dir, "02_expression.rds"))
mat_name <- if (is.null(params$expression)) "tpm" else params$expression
if (is.null(expr[[mat_name]])) {
  stop("expression matrix '", mat_name, "' is not in 02_expression.rds -- re-run the ",
       "expression step, which writes both tpm and tmm")
}
mat  <- expr[[mat_name]]
coh  <- fread(file.path(out_dir, "04_marker_coherence.csv"))
cw   <- fread(file.path(out_dir, "02_sample_crosswalk.csv"))

cts <- if (is.null(params$cell_types)) coh[passes == TRUE, cell_type] else params$cell_types
if (length(cts) < 2L) stop("need at least 2 cell types; got: ", paste(cts, collapse = ", "))

md <- BRETIGEA::markers_df_human_brain
## Genes with zero mean would divide by zero under scaled_mean and contribute
## nothing under the others, so drop them from the pools up front.
nonzero <- rownames(mat)[rowMeans(mat) > 0]
pools <- lapply(cts, function(ct)
  intersect(head(md$markers[md$cell == ct], top_n), nonzero))
names(pools) <- cts
thin <- names(pools)[lengths(pools) < n_mark]
if (length(thin)) {
  stop("cell type(s) with fewer than ", n_mark, " usable markers: ",
       paste(sprintf("%s=%d", thin, lengths(pools)[thin]), collapse = ", "))
}

## ---- stability loop --------------------------------------------------------
## Every aggregator is evaluated on every draw, so the ensemble below is
## available for all of them and not only the primary. Per-gene rescaling does
## not depend on the marker draw, so it is done once here rather than inside
## each of the n_iter calls.
mat_by_agg <- setNames(lapply(aggs, function(a)
                         if (a == "scaled_mean") mat / rowMeans(mat) else mat), aggs)

set.seed(seed)
primary <- aggs[1]
runs <- lapply(seq_len(n_iter), function(i) {
  ml <- lapply(pools, function(g) sample(g, n_mark))
  list(markers = ml,
       prop = setNames(lapply(aggs, function(a)
                dsa_proportions(mat_by_agg[[a]], ml, agg = a, prescaled = TRUE)), aggs))
})
P   <- lapply(runs, function(r) r$prop[[primary]])
n_s <- nrow(P[[1]])

## ---- degenerate cell types -------------------------------------------------
zeroed <- table(unlist(lapply(P, function(p) attr(p, "zeroed"))))
dead   <- names(zeroed)[zeroed == n_iter]
if (length(dead)) {
  stop("cell type(s) zeroed by nnls in every iteration: ", paste(dead, collapse = ", "),
       ".\nTheir marker profiles are not separable from the rest, so DSA cannot ",
       "estimate them here.\nRe-run with --cell-types excluding them (see step 4's ",
       "separation margins).")
}
if (length(zeroed)) {
  cat("\nNOTE: nnls zeroed a cell type in some iterations:\n")
  for (ct in names(zeroed))
    cat(sprintf("  %s: %d of %d iterations\n", ct, zeroed[[ct]], n_iter))
  cat("Those iterations are excluded from that type's stability statistics.\n")
}

## Iterations where a type is degenerate carry no information about it, so they
## are dropped per cell type rather than poisoning the correlation matrix.
stability <- rbindlist(lapply(cts, function(ct) {
  m  <- vapply(P, function(p) p[, ct], numeric(n_s))
  ok <- apply(m, 2, sd) > 0
  if (sum(ok) < 2L) {
    return(data.table(cell_type = ct, median_r = NA_real_, min_r = NA_real_,
                      n_degenerate = sum(!ok)))
  }
  cc   <- suppressWarnings(cor(m[, ok, drop = FALSE]))
  vals <- cc[upper.tri(cc)]
  data.table(cell_type = ct, median_r = median(vals, na.rm = TRUE),
             min_r = min(vals, na.rm = TRUE), n_degenerate = sum(!ok))
}))
stability[, stable := !is.na(median_r) & median_r >= stab_min]

## Representative run: the one whose proportions correlate best, on average
## across cell types, with every other run.
score <- vapply(seq_along(P), function(i) {
  mean(vapply(cts, function(ct) {
    if (sd(P[[i]][, ct]) == 0) return(NA_real_)
    others <- vapply(P[-i], function(q) q[, ct], numeric(n_s))
    ok <- apply(others, 2, sd) > 0
    if (!any(ok)) return(NA_real_)
    mean(suppressWarnings(cor(P[[i]][, ct], others[, ok, drop = FALSE])))
  }, numeric(1)), na.rm = TRUE)
}, numeric(1))
if (all(is.na(score))) stop("no iteration produced a usable estimate for any cell type")
rep_idx <- which.max(replace(score, is.na(score), -Inf))
final_markers <- runs[[rep_idx]]$markers

## ---- ensemble estimate -----------------------------------------------------
## The n_iter draws exist to establish that the estimate does not depend on
## which markers were sampled; the between-run Pearson correlations above are
## that check. Once it passes, the ensemble mean is a strictly better point
## estimate than any single draw, and the between-iteration SD is the only
## measure of marker-sampling uncertainty this design produces. Steps 6 and 8
## treat these proportions as known covariates, so that SD needs to be on the
## record -- error in a regressor attenuates its coefficient, which in step 8
## means understating exactly the confounding the analysis is there to measure.
##
## Iterations where nnls zeroed a cell type carry no information about it, so
## they are dropped per type, the same rule the stability table uses. Each
## draw's rows already sum to 1, so the mean does too and the renormalisation
## below is a no-op except where those drops were uneven.
ensemble_of <- function(agg_name) {
  Pa <- lapply(runs, function(r) r$prop[[agg_name]])
  mu <- sdv <- matrix(NA_real_, n_s, length(cts),
                      dimnames = list(rownames(Pa[[1]]), cts))
  for (ct in cts) {
    v  <- vapply(Pa, function(p) p[, ct], numeric(n_s))
    ok <- apply(v, 2, sd) > 0
    if (!any(ok)) next
    mu[, ct]  <- rowMeans(v[, ok, drop = FALSE])
    sdv[, ct] <- apply(v[, ok, drop = FALSE], 1, sd)
  }
  ## na.rm guards the case where a type was never estimable and its column
  ## stayed NA; the `dead` check above should already have stopped, but a whole
  ## proportion matrix silently becoming NA is not a failure worth risking.
  list(mean = mu / rowSums(mu, na.rm = TRUE), sd = sdv)
}

## ---- which estimate is primary ---------------------------------------------
## Averaging is conditional on the stability gate already computed above: if
## the marker draw demonstrably moves a cell type, averaging over draws hides
## that rather than resolving it.
##
## Named prop_mode and est_used rather than mode and estimator: `mode` is a
## base function, and `estimator` is also a column name below, where a bare
## RHS reference inside `:=` would depend on the column not yet existing.
prop_mode <- if (is.null(params$proportions)) "auto" else params$proportions
if (!prop_mode %in% c("auto", "ensemble", "representative"))
  stop("--proportions must be one of: auto, ensemble, representative")
unstable <- stability[stable == FALSE, cell_type]
use_ensemble <- switch(prop_mode, ensemble = TRUE, representative = FALSE,
                       auto = length(unstable) == 0L)
est_used <- if (use_ensemble) "ensemble" else "representative"

if (prop_mode == "auto" && !use_ensemble) {
  cat(sprintf("\nWARNING: %s did not clear median between-run r >= %.2f, so proportions\n",
              paste(unstable, collapse = ", "), stab_min))
  cat("come from the representative run rather than the ensemble mean. Drop the type\n")
  cat("with --cell-types and re-run, or force the ensemble with --proportions ensemble\n")
  cat("and report the instability alongside it.\n")
}

## ---- final proportions, all aggregators ------------------------------------
long <- rbindlist(lapply(aggs, function(a) {
  if (use_ensemble) {
    e <- ensemble_of(a); p <- e$mean; s <- e$sd
  } else {
    p <- dsa_proportions(mat_by_agg[[a]], final_markers, agg = a, prescaled = TRUE)
    s <- matrix(NA_real_, nrow(p), ncol(p), dimnames = dimnames(p))
  }
  m1 <- melt(as.data.table(p, keep.rownames = "rids"), id.vars = "rids",
             variable.name = "cell_type", value.name = "proportion")
  m2 <- melt(as.data.table(s, keep.rownames = "rids"), id.vars = "rids",
             variable.name = "cell_type", value.name = "proportion_sd")
  merge(m1, m2, by = c("rids", "cell_type"))[, `:=`(aggregator = a,
                                                    estimator = est_used)]
}))
long <- merge(long, cw[, .(rids, con)], by = "rids")
setcolorder(long, c("rids", "con", "cell_type", "aggregator", "estimator",
                    "proportion", "proportion_sd"))
fwrite(long, outf("cell_proportions", "csv"))

fwrite(stability, outf("dsa_stability", "csv"))
## The representative run's markers, kept as a record. Under the ensemble no
## single marker set defines the estimate, so this is provenance rather than
## the definition of what was computed.
fwrite(rbindlist(lapply(names(final_markers), function(ct)
         data.table(cell_type = ct, marker = final_markers[[ct]]))),
       outf("dsa_markers", "csv"))

## ---- plot ------------------------------------------------------------------
prim <- long[aggregator == primary]
png(outf("proportions", "png"), width = 1400, height = 900, res = 150)
boxplot(proportion ~ con + cell_type, data = prim, las = 2, cex.axis = 0.7,
        col = c("grey80", "steelblue"), ylab = "estimated relative proportion",
        xlab = "", main = sprintf("DSA proportions (%s aggregation)", primary))
invisible(dev.off())

## ---- report ----------------------------------------------------------------
cat("\n== DSA ==\n")
cat(sprintf("expression matrix  : %s\n", mat_name))
cat(sprintf("cell types         : %s\n", paste(cts, collapse = ", ")))
cat(sprintf("markers per type   : %d sampled from pools of %s\n",
            n_mark, paste(lengths(pools), collapse = "/")))
cat(sprintf("iterations / seed  : %d / %d\n", n_iter, seed))
cat(sprintf("representative run : %d\n", rep_idx))
cat(sprintf("estimator          : %s%s\n", est_used,
            if (use_ensemble) sprintf(" (mean of %d draws)", n_iter)
            else sprintf(" (run %d)", rep_idx)))

cat("\n== Between-run stability ==\n")
print(stability[, .(cell_type, median_r = round(median_r, 3),
                    min_r = round(min_r, 3), n_degenerate, stable)], row.names = FALSE)
if (any(!stability$stable)) {
  cat(sprintf("\nCAUTION: median between-run r < %.2f for %s -- report with explicit\n",
              stab_min, paste(stability[stable == FALSE, cell_type], collapse = ", ")))
  cat("caution or drop, as Wang et al. did for OPC.\n")
}

if (use_ensemble) {
  cat("\n== Between-iteration SD of the ensemble (across samples) ==\n")
  print(long[aggregator == primary,
             .(median_sd = round(median(proportion_sd, na.rm = TRUE), 4),
               max_sd    = round(max(proportion_sd, na.rm = TRUE), 4)),
             by = cell_type], row.names = FALSE)
  cat("This is marker-sampling uncertainty only. It is NOT the uncertainty in\n")
  cat("these proportions as estimates of true composition, which is larger.\n")
}

cat(sprintf("\n== Median proportion by group (%s) ==\n", primary))
print(dcast(prim, cell_type ~ con, value.var = "proportion",
            fun.aggregate = function(x) round(median(x), 4)), row.names = FALSE)

write_json(list(expression = mat_name,
                cell_types = cts, top_n = top_n, n_markers = n_mark,
                n_iter = n_iter, seed = seed, rep_iteration = rep_idx,
                aggregators = aggs, primary_aggregator = primary,
                estimator = est_used, proportions_mode = prop_mode,
                unstable_cell_types = unstable,
                stability_min = stab_min, stability = stability,
                label = params$label,
                n_zeroed_by_nnls = as.list(zeroed),
                marker_pool_sizes = as.list(lengths(pools))),
           outf("dsa_summary", "json"), auto_unbox = TRUE, pretty = TRUE)

cat("\nProportions are RELATIVE shares of the tissue mRNA pool, not cell counts.\n")
cat(sprintf("Wrote results to %s\n", out_dir))
