#!/usr/bin/env Rscript
##
## 06_stats.R -- group comparison, neuropathology correlation and covariate
## adjustment for the estimated proportions.
##
##   Rscript 06_stats.R <params.json>
##
## Confidence intervals are the point. At modest n a non-significant result is
## "consistent with a small change in either direction", not "no change".
##

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(betareg); library(effsize)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 06_stats.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir   <- params$output_dir
sfx       <- if (is.null(params$label) || !nzchar(params$label)) "" else paste0("_", params$label)
cont_vars <- params$continuous
cov_vars  <- params$covariates
med_vars  <- if (is.null(params$mediators)) character(0) else params$mediators
primary   <- params$aggregator
outf <- function(stem, ext) file.path(out_dir, sprintf("06_%s%s.%s", stem, sfx, ext))

props <- fread(file.path(out_dir, sprintf("05_cell_proportions%s.csv", sfx)))
covs  <- fread(file.path(out_dir, "06_covariates.csv"))
d <- merge(props, covs[, setdiff(names(covs), "con"), with = FALSE], by = "rids")

groups <- sort(unique(d$con))
if (length(groups) != 2L) stop("expected exactly 2 groups in `con`, found: ",
                               paste(groups, collapse = ", "))
case <- if (is.null(params$case_group)) groups[1] else params$case_group
if (!case %in% groups) stop("--case-group ", case, " is not one of: ",
                            paste(groups, collapse = ", "))
ctrl <- setdiff(groups, case)
cat(sprintf("\nTreating '%s' as case and '%s' as reference.\n", case, ctrl))
cat("Every ratio, percent change and shift below is case relative to reference.\n")

## ---- binary comparison -----------------------------------------------------
## A degenerate cell type -- nnls pinned to zero in every sample, which happens
## when a type is not separable under a given aggregator -- has no variance for
## the rank test to work on, and wilcox.test's Hodges-Lehmann root-finder fails
## outright rather than returning NA. The p-value is still defined, so only the
## interval is dropped: the cell type stays in the table marked NA instead of
## aborting the step for every other type.
hl_test <- function(x, y) {
  wt <- tryCatch(suppressWarnings(wilcox.test(x, y, conf.int = TRUE, exact = FALSE)),
                 error = function(e) NULL)
  if (!is.null(wt) && !is.null(wt$conf.int)) {
    return(list(hl = unname(wt$estimate), lo = wt$conf.int[1],
                hi = wt$conf.int[2], p = wt$p.value))
  }
  wt <- tryCatch(suppressWarnings(wilcox.test(x, y, exact = FALSE)),
                 error = function(e) NULL)
  list(hl = NA_real_, lo = NA_real_, hi = NA_real_,
       p = if (is.null(wt)) NA_real_ else wt$p.value)
}

binary <- d[, {
  x <- proportion[con == case]; y <- proportion[con == ctrl]
  wt <- hl_test(x, y)
  cd <- tryCatch(suppressWarnings(effsize::cliff.delta(x, y)),
                 error = function(e) list(estimate = NA_real_))
  mr <- median(y)
  .(n_case = length(x), n_ref = length(y),
    median_case = median(x), median_ref = mr,
    q1_case = quantile(x, .25), q3_case = quantile(x, .75),
    q1_ref  = quantile(y, .25), q3_ref  = quantile(y, .75),
    ## A zero reference median leaves the ratio undefined rather than infinite:
    ## a type that could not be estimated has no fold change to report.
    ratio = if (mr > 0) median(x) / mr else NA_real_,
    pct_change = if (mr > 0) 100 * (median(x) - mr) / mr else NA_real_,
    hl = wt$hl, hl_lo = wt$lo, hl_hi = wt$hi,
    pval = wt$p, cliffs_delta = unname(cd$estimate))
}, by = .(cell_type, aggregator)]
binary[, fdr := p.adjust(pval, "BH"), by = aggregator]
fwrite(binary, outf("stats_binary", "csv"))

## ---- continuous: does composition track pathology? -------------------------
## Spearman CI via the Fisher z transform with the standard 1.06 inflation for
## rank correlations. Approximate, but adequate beside a rank-based test.
cont <- rbindlist(lapply(cont_vars, function(v) {
  if (!v %in% names(d)) {
    warning("continuous variable not found, skipping: ", v); return(NULL)
  }
  d[, {
    yy <- suppressWarnings(as.numeric(get(v)))
    ok <- is.finite(proportion) & is.finite(yy)
    if (sum(ok) < 5L || length(unique(yy[ok])) < 3L) {
      .(variable = v, n = sum(ok), rho = NA_real_, lo = NA_real_,
        hi = NA_real_, pval = NA_real_)
    } else {
      ct <- suppressWarnings(cor.test(proportion[ok], yy[ok], method = "spearman"))
      n  <- sum(ok); z <- atanh(unname(ct$estimate)); se <- 1.06 / sqrt(n - 3)
      .(variable = v, n = n, rho = unname(ct$estimate),
        lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se), pval = ct$p.value)
    }
  }, by = .(cell_type, aggregator)]
}))
if (nrow(cont)) {
  cont[, fdr := p.adjust(pval, "BH"), by = .(aggregator, variable)]
  fwrite(cont, outf("stats_continuous", "csv"))
}

## ---- adjusted model --------------------------------------------------------
## Proportions are bounded in (0,1), so beta regression rather than OLS; values
## are squeezed off the boundary because betareg excludes exact 0 and 1.
d[, diagnosis := factor(con, levels = c(ctrl, case))]
## betareg excludes exact 0 and 1, which nnls can produce. Smithson and
## Verkuilen's (y*(n-1) + 0.5)/n compresses toward 0.5 by an amount that scales
## with the sample size, rather than pinning boundary values to a constant
## chosen by hand: with an arbitrary epsilon, a cell type carrying several exact
## zeros has its coefficient determined by the value of that epsilon.
d[, y := (proportion * (.N - 1) + 0.5) / .N, by = .(cell_type, aggregator)]

## Some cohorts censor age at death as a string such as "90+". Left as-is it
## becomes NA on coercion and those participants drop silently out of every
## adjusted model.
if ("age_death" %in% names(d) && !is.numeric(d$age_death)) {
  raw <- as.character(d$age_death)
  n_censored <- sum(grepl("\\+", raw), na.rm = TRUE)
  d[, age_death := suppressWarnings(as.numeric(sub("\\+.*$", "", raw)))]
  if (n_censored > 0L) {
    cat(sprintf("age_death: %d row(s) censored as \"90+\" recoded to the numeric bound.\n",
                n_censored))
  }
  if (any(is.na(d$age_death))) {
    cat(sprintf("WARNING: age_death is NA for %d row(s) after coercion.\n",
                sum(is.na(d$age_death))))
  }
}

## A covariate perfectly collinear with the grouping leaves the model
## unidentifiable -- it has no within-group variation to contribute.
screen_covariates <- function(vars, label) {
  keep <- character(0)
  for (v in vars) {
    if (!v %in% names(d)) { cat(sprintf("  %s: not found, skipping\n", v)); next }
    vals <- d[[v]]
    if (length(unique(vals[!is.na(vals)])) < 2L) {
      cat(sprintf("  %s: constant, skipping\n", v)); next
    }
    tab <- table(d$con, vals)
    if (all(rowSums(tab > 0) == 1L) && ncol(tab) >= 2L) {
      cat(sprintf("  %s: perfectly collinear with `con`, dropped from the %s model\n",
                  v, label))
      next
    }
    keep <- c(keep, v)
  }
  keep
}

fit_adjusted <- function(vars, model_name) {
  form <- as.formula(paste("y ~ diagnosis",
                           if (length(vars)) paste("+", paste(vars, collapse = " + ")) else ""))
  out <- rbindlist(lapply(split(d, by = c("cell_type", "aggregator")), function(sub) {
    fit <- tryCatch(betareg(form, data = sub), error = function(e) NULL)
    coefname <- paste0("diagnosis", case)
    if (is.null(fit) || !coefname %in% rownames(coef(summary(fit))$mean)) {
      return(data.table(model = model_name, cell_type = sub$cell_type[1],
                        aggregator = sub$aggregator[1], estimate = NA_real_,
                        std_error = NA_real_, pval = NA_real_, converged = FALSE))
    }
    cf <- coef(summary(fit))$mean[coefname, ]
    data.table(model = model_name, cell_type = sub$cell_type[1],
               aggregator = sub$aggregator[1], estimate = cf[["Estimate"]],
               std_error = cf[["Std. Error"]], pval = cf[["Pr(>|z|)"]], converged = TRUE)
  }), fill = TRUE)
  out[converged == TRUE, fdr := p.adjust(pval, "BH"), by = aggregator]
  out[]
}

## ---- covariate screen ------------------------------------------------------
## A covariate earns its parameter only if it plausibly confounds: it must both
## differ between groups AND predict the outcome. At n = 59 each predictor costs
## real precision, so the screen is reported rather than left implicit. It is
## NOT applied automatically -- data-dependent model selection has its own
## problems, and the decision belongs to the analyst via --covariates.
wide_props <- dcast(d[aggregator == primary], rids + con ~ cell_type,
                    value.var = "proportion")
scr_input <- merge(wide_props,
                   unique(d[, c("rids", intersect(c(cov_vars, med_vars), names(d))),
                            with = FALSE]), by = "rids")
cts_here <- setdiff(names(wide_props), c("rids", "con"))

screen <- rbindlist(lapply(intersect(c(cov_vars, med_vars), names(scr_input)), function(v) {
  vals <- suppressWarnings(as.numeric(scr_input[[v]]))
  if (all(is.na(vals)) || length(unique(vals[!is.na(vals)])) < 2L) return(NULL)
  ## A rank test on a binary covariate is not the right instrument -- msex and
  ## the like are contingency tables, so they get Fisher's exact.
  uv <- unique(vals[!is.na(vals)])
  gw_p <- if (length(uv) <= 2L) {
    tryCatch(fisher.test(table(scr_input$con, vals))$p.value,
             error = function(e) NA_real_)
  } else {
    tryCatch(suppressWarnings(wilcox.test(vals ~ scr_input$con, exact = FALSE))$p.value,
             error = function(e) NA_real_)
  }
  rows <- data.table(covariate = v, target = "diagnosis",
                     statistic = NA_real_, pval = gw_p)
  rows <- rbind(rows, rbindlist(lapply(cts_here, function(ct) {
    ct_test <- suppressWarnings(cor.test(vals, scr_input[[ct]], method = "spearman"))
    data.table(covariate = v, target = ct,
               statistic = unname(ct_test$estimate), pval = ct_test$p.value)
  })))
  rows[]
}))
if (nrow(screen)) {
  fwrite(screen, outf("covariate_screen", "csv"))
  cat("\n== Covariate screen ==\n")
  print(dcast(screen, covariate ~ target, value.var = "pval",
              fun.aggregate = function(x) signif(min(x), 3)), row.names = FALSE)
  cat("p-values: association with diagnosis, and Spearman against each proportion.\n")
  cat("A covariate null on both is not earning its parameter; null on diagnosis but\n")
  cat("associated with composition still buys precision. Reported, not applied --\n")
  cat("adjust --covariates yourself if the screen argues for it.\n")
}

cat("\nScreening primary covariates:\n")
usable <- screen_covariates(cov_vars, "primary")
cat(sprintf("Primary adjusted model: y ~ diagnosis%s\n",
            if (length(usable)) paste0(" + ", paste(usable, collapse = " + ")) else ""))
adjusted <- fit_adjusted(usable, "primary")

## Mediator sensitivity model, reported separately and never as "the adjusted
## result". A variable on the causal path from diagnosis to composition removes
## part of the very effect being estimated when conditioned on.
med_usable <- character(0)
if (length(med_vars)) {
  cat("\nScreening mediator covariates (sensitivity model only):\n")
  med_usable <- screen_covariates(med_vars, "mediator")
  if (length(med_usable)) {
    cat(sprintf("Mediator sensitivity model: y ~ diagnosis%s + %s\n",
                if (length(usable)) paste0(" + ", paste(usable, collapse = " + ")) else "",
                paste(med_usable, collapse = " + ")))
    adjusted <- rbind(adjusted, fit_adjusted(c(usable, med_usable), "mediator_sensitivity"),
                      fill = TRUE)
  }
}
fwrite(adjusted, outf("stats_adjusted", "csv"))

## ---- forest plot -----------------------------------------------------------
## The central point of the whole analysis: each interval is a statement about
## how much change these data are consistent with, in either direction. A table
## states that; a forest plot shows it.
##
## This figure carries THIS cohort's estimates only. Effects from other
## cohorts are deliberately not drawn on it: a reference marker on the same
## axis reads as a comparison this analysis did not make, and invites the
## interval to be judged against a cross-cohort, cross-method number it was
## never tested against. Any such comparison belongs in prose, where its
## provenance can be stated.
fp <- binary[aggregator == primary][order(cell_type)]

png(outf("forest", "png"), width = 1200, height = 800, res = 150)
par(mar = c(5, 6, 4, 2))
yy <- seq_len(nrow(fp))
xlim <- range(c(fp$hl_lo, fp$hl_hi, 0), na.rm = TRUE)
plot(fp$hl, yy, xlim = xlim, ylim = c(0.5, nrow(fp) + 0.5), yaxt = "n",
     pch = 19, cex = 1.3, xlab = "Hodges-Lehmann shift in proportion (case - reference)",
     ylab = "", main = sprintf("Cell-type composition: %s vs %s", case, ctrl))
axis(2, at = yy, labels = fp$cell_type, las = 1)
segments(fp$hl_lo, yy, fp$hl_hi, yy, lwd = 2)
abline(v = 0, lty = 1, col = "grey40")
legend("topright", "observed (95% CI)", pch = 19, col = "black",
       bty = "n", cex = 0.8)
invisible(dev.off())

## ---- plot ------------------------------------------------------------------
prim <- d[aggregator == primary]
## Display labels for the two groups. The database values are used everywhere
## else in this script; this is presentation only.
glab     <- setNames(c("AD", "CN"), c(case, ctrl))
lev      <- levels(factor(prim$con))
cts_plot <- sort(unique(as.character(prim$cell_type)))
ytop     <- max(prim$proportion, na.rm = TRUE)

png(outf("proportions_by_group", "png"), width = 1400, height = 900, res = 150)
par(mar = c(5, 5, 2, 2))
boxplot(proportion ~ con + cell_type, data = prim, outline = FALSE, xaxt = "n",
        col = c("grey80", "steelblue"),
        ylab = "Estimated relative proportion", xlab = "Cell type",
        ylim = c(0, ytop * 1.15))
stripchart(proportion ~ con + cell_type, data = prim, vertical = TRUE, add = TRUE,
           method = "jitter", jitter = 0.15, pch = 19, cex = 0.5,
           col = adjustcolor("black", 0.4))
## One tick per cell type, centred between that type's pair of boxes. Group is
## carried by the legend instead, so the axis is not repeating it ten times.
at_ct <- seq(1.5, by = 2, length.out = length(cts_plot))
axis(1, at = at_ct, labels = cts_plot, tick = FALSE)
legend("topright", legend = unname(glab[lev]), fill = c("grey80", "steelblue"),
       bty = "n", cex = 0.85)
invisible(dev.off())

## ---- report ----------------------------------------------------------------
cat(sprintf("\n== Binary comparison (%s) ==\n", primary))
print(binary[aggregator == primary,
             .(cell_type, median_case = round(median_case, 4),
               median_ref = round(median_ref, 4), ratio = round(ratio, 3),
               pct = round(pct_change, 1), pval = signif(pval, 3),
               fdr = signif(fdr, 3))], row.names = FALSE)

cat(sprintf("\n== Effect size with CI (%s) ==\n", primary))
print(binary[aggregator == primary,
             .(cell_type, hl = signif(hl, 3),
               CI = paste0(signif(hl_lo, 3), " to ", signif(hl_hi, 3)),
               cliffs_delta = round(cliffs_delta, 3))], row.names = FALSE)
cat("A CI spanning zero means the data are consistent with a change in either\n")
cat("direction -- report it that way, not as an absence of change.\n")

if (nrow(cont)) {
  cat(sprintf("\n== Correlation with neuropathology (%s) ==\n", primary))
  print(cont[aggregator == primary,
             .(cell_type, variable, n, rho = round(rho, 3),
               CI = paste0(round(lo, 2), " to ", round(hi, 2)),
               pval = signif(pval, 3), fdr = signif(fdr, 3))], row.names = FALSE)
}

cat(sprintf("\n== Adjusted diagnosis effect, PRIMARY model (%s, logit scale) ==\n", primary))
print(adjusted[model == "primary" & aggregator == primary,
               .(cell_type, estimate = round(estimate, 3),
                 std_error = round(std_error, 3), pval = signif(pval, 3),
                 fdr = signif(fdr, 3), converged)], row.names = FALSE)

if (any(adjusted$model == "mediator_sensitivity")) {
  cat(sprintf("\n== Mediator sensitivity model (%s) -- NOT the adjusted result ==\n", primary))
  print(adjusted[model == "mediator_sensitivity" & aggregator == primary,
                 .(cell_type, estimate = round(estimate, 3),
                   std_error = round(std_error, 3), pval = signif(pval, 3),
                   fdr = signif(fdr, 3))], row.names = FALSE)
  cat(sprintf("Adds %s, which is plausibly on the causal path from diagnosis to\n",
              paste(med_usable, collapse = " + ")))
  cat("composition. Conditioning on it removes part of the effect being measured,\n")
  cat("so a coefficient that shrinks or changes sign here is expected behaviour,\n")
  cat("not evidence about the diagnosis effect. Report alongside the primary\n")
  cat("model with that caveat, never in place of it.\n")
}

write_json(list(label = params$label, case_group = case, reference_group = ctrl,
                primary_aggregator = primary, continuous = cont_vars,
                covariates_used = usable, covariates_requested = cov_vars,
                mediator_covariates_used = med_usable,
                mediator_covariates_requested = med_vars,
                binary = binary, continuous_results = cont, adjusted = adjusted),
           outf("stats_summary", "json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("\nWrote results to %s\n", out_dir))
