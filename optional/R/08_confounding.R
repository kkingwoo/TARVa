#!/usr/bin/env Rscript
##
## 08_confounding.R -- does estimated cell composition account for the calls
## that separate the two groups?
##
##   Rscript 08_confounding.R <params.json>
##
## THE MECHANISM. Calls are made only where dp >= dp_threshold. So a composition
## shift changes expression, changes depth, changes how many positions clear the
## threshold, and can make a gene look condition-enriched with no change in
## editing biology at all. Presence/absence findings are the most exposed.
##
## WHAT THIS DATA CAN AND CANNOT SETTLE. dp is recorded only where a call was
## made. A (gene, sample) pair with no call is therefore ambiguous: it is either
## a real negative, or a position that never had the depth to be tested. Nothing
## here can tell those apart -- that needs the VCFs or BAMs, which no longer
## exist for this cohort.
##
## The response is to adjust for the OBSERVED determinants of depth rather than
## to impute which zeros were testable. Gene expression is observed for every
## (gene, sample) pair, and library size for every sample; both are carried as
## covariates. Because every model below is fitted WITHIN a gene, anything
## affecting coverage that is gene-specific but sample-constant -- 3' bias,
## transcript length, mappability -- is absorbed by the fit. What escapes is
## sample-by-gene interaction in coverage bias, which is second-order and has no
## reason to align with diagnosis. That is the residual, and it is named rather
## than hidden.
##
## Nothing in this script extrapolates depth to positions where it was not
## observed. An earlier version scaled one group's median depth by an expression
## ratio to argue a position "was testable"; that fits on called positions and
## applies the result to uncalled ones, which is the selection problem this
## analysis exists to investigate. It is gone.
##

suppressPackageStartupMessages({
  library(jsonlite); library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 08_confounding.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir  <- params$output_dir
sfx      <- if (is.null(params$label) || !nzchar(params$label)) "" else paste0("_", params$label)
agg      <- params$aggregator
seed     <- params$seed
min_n        <- if (is.null(params$min_n))        20L   else params$min_n
min_carriers <- if (is.null(params$min_carriers)) 5L    else params$min_carriers
dp_thresh    <- if (is.null(params$dp_threshold)) 20L   else params$dp_threshold
r2_margin    <- if (is.null(params$r2_margin))    0.10  else params$r2_margin
n_boot       <- if (is.null(params$n_bootstrap))  2000L else params$n_bootstrap
outf <- function(stem, ext) file.path(out_dir, sprintf("08_%s%s.%s", stem, sfx, ext))

props  <- fread(file.path(out_dir, sprintf("05_cell_proportions%s.csv", sfx)))[aggregator == agg]
burden <- fread(file.path(out_dir, "08_burden.csv"))
sites  <- fread(file.path(out_dir, "08_gene_sites.csv"))
groups <- fread(file.path(out_dir, "08_gene_group_resolved.csv"))
expr   <- readRDS(file.path(out_dir, "02_expression.rds"))

## ---- which expression matrix -----------------------------------------------
## Match the normalisation the deconvolution actually ran on, the way step 7
## does. Modelling composition against a differently-normalised matrix than the
## one the proportions came from puts a per-sample scale difference into every
## coefficient. The *_all matrices are used because they carry the labelled
## genes filterByExpr dropped -- disproportionately the low-expressed ones
## nearest the detection floor, which are exactly the genes most exposed to the
## mechanism. Step 2 builds them on the filtered per-sample denominators, so
## they are a strict superset rather than a second scale.
dsa_sum  <- file.path(out_dir, sprintf("05_dsa_summary%s.json", sfx))
mat_name <- if (!is.null(params$expression)) params$expression else
            if (file.exists(dsa_sum)) fromJSON(dsa_sum)$expression else "tpm"
all_name <- paste0(mat_name, "_all")
if (is.null(expr[[all_name]])) {
  stop("expression matrix '", all_name, "' is not in 02_expression.rds -- re-run the ",
       "expression step, which writes the unfiltered matrices step 8 needs")
}
mat <- expr[[all_name]]
cat(sprintf("Modelling against the %s matrix (%d genes, unfiltered).\n",
            all_name, nrow(mat)))

## Library size, an observed determinant of depth. Counts are the filtered
## matrix, which is the right denominator: it is the sequencing effort that
## reached analysable genes.
lib_size <- colSums(expr$counts)

wide <- dcast(props, rids + con ~ cell_type, value.var = "proportion")
cts  <- setdiff(names(wide), c("rids", "con"))

## Every model below joins expression, proportions and library size on rids. A
## silent mismatch between the matrix and the proportions file would not error,
## it would quietly shrink the sample -- so it is checked rather than assumed,
## and the matrix is put in the proportions file's order once for all of them.
miss_prop <- setdiff(colnames(mat), wide$rids)
miss_mat  <- setdiff(wide$rids, colnames(mat))
if (length(miss_prop) || length(miss_mat)) {
  stop("sample mismatch between 02_expression.rds and the proportions file.\n",
       if (length(miss_prop))
         paste0("  in the matrix but not the proportions: ",
                paste(head(miss_prop, 10), collapse = ", "), "\n") else "",
       if (length(miss_mat))
         paste0("  in the proportions but not the matrix: ",
                paste(head(miss_mat, 10), collapse = ", "), "\n") else "",
       "Re-run the expression and dsa steps into the same --output.")
}
mat      <- mat[, wide$rids, drop = FALSE]
lib_size <- lib_size[wide$rids]

grps <- sort(unique(wide$con))
case <- if (is.null(params$case_group)) grps[1] else params$case_group
ctrl <- setdiff(grps, case)
cat(sprintf("Treating '%s' as case, '%s' as reference.\n", case, ctrl))

## Proportions sum to 1, so one cell type must come out of every model.
drop_ct <- if (!is.null(params$drop_cell_type)) params$drop_cell_type else
             cts[which.max(vapply(cts, function(c) mean(wide[[c]]), numeric(1)))]
keep_ct <- setdiff(cts, drop_ct)
cat(sprintf("Dropping '%s' from all models (most abundant); modelling on: %s\n",
            drop_ct, paste(keep_ct, collapse = ", ")))

## ---- shared lookups --------------------------------------------------------
## Built once here because 8b-norm, 8f and 8g all need them.
labels <- unique(rbind(groups[, .(ref_base = as.character(gene),   gene_group)],
                       groups[, .(ref_base = as.character(ensgid), gene_group)]))
labels <- labels[!is.na(ref_base) & nzchar(ref_base)]
## One row per key. Without this a key appearing under two class labels -- or a
## gene and ensgid colliding -- duplicates rows on every merge below, silently
## inflating counts and the percentages computed from them.
n_dup_keys <- sum(duplicated(labels$ref_base))
if (n_dup_keys) {
  cat(sprintf("NOTE: %d gene key(s) carry more than one class label; keeping the first.\n",
              n_dup_keys))
  labels <- labels[!duplicated(ref_base)]
}
sym_of <- c(setNames(groups$resolved_symbol, as.character(groups$gene)),
            setNames(groups$resolved_symbol, as.character(groups$ensgid)))

## Long-form expression for the labelled genes, one row per (gene, sample).
## This is the object that makes the ambiguous zeros tractable: it exists for
## every pair, whether or not a call was made.
gene_keys <- unique(sites$ref_base)
expr_long <- rbindlist(lapply(gene_keys, function(k) {
  s <- unname(sym_of[as.character(k)])
  if (is.na(s) || !s %in% rownames(mat)) return(NULL)
  data.table(ref_base = k, rids = colnames(mat), expression = as.numeric(mat[s, ]))
}))
cat(sprintf("Expression resolved for %d of %d gene keys carrying calls.\n",
            uniqueN(expr_long$ref_base), length(gene_keys)))

## ============================================================================
## 8a -- per-individual call burden
## ============================================================================
## Promoted to the front because it is the only component with no missing-data
## problem at all. Burden is a property of the individual, counted over
## everything, and is not conditional on any gene being testable. If composition
## does not predict how many calls a person yields, the mechanism has no
## traction at scale -- a stronger result than any per-gene adjustment, because
## it addresses the cause rather than the traces.
b <- merge(burden, wide, by = "rids", suffixes = c("", "_p"))
bx <- b[con == case, n_calls]; by_ <- b[con == ctrl, n_calls]
bw <- suppressWarnings(wilcox.test(bx, by_, conf.int = TRUE, exact = FALSE))
per_capita <- c(case = mean(bx), ref = mean(by_))

## Marginal AND partial. Composition and burden both track diagnosis, so the
## marginal correlation conflates the two; the partial is Spearman on ranks
## residualised against diagnosis, which needs no extra package and no
## distributional assumption.
partial_spearman <- function(x, y, g) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L) return(c(NA_real_, NA_real_))
  rx <- resid(lm(rank(x[ok]) ~ factor(g[ok])))
  ry <- resid(lm(rank(y[ok]) ~ factor(g[ok])))
  ct <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
  c(unname(ct$estimate), ct$p.value)
}

burden_cmp <- rbindlist(lapply(cts, function(ct) {
  m <- suppressWarnings(cor.test(b$n_calls, b[[ct]], method = "spearman"))
  p <- partial_spearman(b$n_calls, b[[ct]], b$con)
  data.table(cell_type = ct, rho = unname(m$estimate), pval = m$p.value,
             rho_partial = p[1], pval_partial = p[2])
}))
burden_cmp[, `:=`(fdr = p.adjust(pval, "BH"), fdr_partial = p.adjust(pval_partial, "BH"))]
fwrite(burden_cmp, outf("burden_vs_composition", "csv"))

cat("\n== 8a. Per-individual call burden (no missing-data problem) ==\n")
cat(sprintf("mean calls per individual: %s %.1f, %s %.1f  (ratio %.2f)\n",
            case, per_capita[["case"]], ctrl, per_capita[["ref"]],
            per_capita[["case"]] / per_capita[["ref"]]))
cat(sprintf("group difference: p = %.3g\n\n", bw$p.value))
print(burden_cmp[, .(cell_type, rho = round(rho, 3), fdr = signif(fdr, 3),
                     rho_partial = round(rho_partial, 3),
                     fdr_partial = signif(fdr_partial, 3))], row.names = FALSE)
cat("rho is marginal; rho_partial removes diagnosis from both sides. Burden\n")
cat("unrelated to composition means the mechanism has no foothold.\n")

## ============================================================================
## 8b -- per-gene: how much of expression does composition explain?
## ============================================================================
groups_m <- groups[!is.na(resolved_symbol) & resolved_symbol %in% rownames(mat)]
cat(sprintf("\n%d labelled genes have expression and are modelled.\n", nrow(groups_m)))

pm <- as.matrix(wide[, keep_ct, with = FALSE])
rownames(pm) <- wide$rids
pm <- pm[colnames(mat), , drop = FALSE]

gene_models <- rbindlist(lapply(seq_len(nrow(groups_m)), function(i) {
  sym <- groups_m$resolved_symbol[i]
  y <- log2(mat[sym, ] + 1)
  if (sd(y) == 0) return(NULL)
  fit <- lm(y ~ pm)
  cf  <- coef(fit)
  out <- data.table(gene = groups_m$gene[i], gene_name = sym,
                    gene_group = groups_m$gene_group[i],
                    r2 = summary(fit)$r.squared,
                    ## Adjusted R2 is the one to compare and to read against the
                    ## equivalence margin: with k predictors on n samples the raw
                    ## statistic has an expectation of about k/(n-1) under the
                    ## null, which is not far below the margin itself.
                    adj_r2 = summary(fit)$adj.r.squared,
                    pval = tryCatch({
                      f <- summary(fit)$fstatistic
                      unname(pf(f[1], f[2], f[3], lower.tail = FALSE))
                    }, error = function(e) NA_real_))
  for (ct in keep_ct) out[[paste0("beta_", ct)]] <- unname(cf[paste0("pm", ct)])
  out
}))
gene_models[, fdr := p.adjust(pval, "BH")]
fwrite(gene_models, outf("gene_models", "csv"))

null_r2 <- length(keep_ct) / (nrow(wide) - 1)
cat("\n== 8b. Composition R2 per gene, by class ==\n")
print(gene_models[, .(n = .N,
                      median_r2     = round(median(r2), 3),
                      median_adj_r2 = round(median(adj_r2), 3),
                      q1 = round(quantile(adj_r2, .25), 3),
                      q3 = round(quantile(adj_r2, .75), 3),
                      pct_fdr05 = round(100 * mean(fdr < 0.05), 1)),
                  by = gene_group], row.names = FALSE)
cat(sprintf("Raw R2 has an expectation of about %.3f under the null here (%d predictors,\n",
            null_r2, length(keep_ct)))
cat(sprintf("%d samples); adj_r2 removes that floor and is the column to read.\n", nrow(wide)))
cat("This is how much of these genes' EXPRESSION composition explains. It is a\n")
cat("proxy for confounding of the calls, not a measurement of it -- 8f is that.\n")
kw <- suppressWarnings(kruskal.test(adj_r2 ~ factor(gene_group), data = gene_models))
cat(sprintf("Kruskal-Wallis across classes: p = %.3g\n", kw$p.value))

## ---- equivalence, not just non-significance ---------------------------------
set.seed(seed)
classes <- sort(unique(gene_models$gene_group))
equiv <- rbindlist(lapply(combn(classes, 2, simplify = FALSE), function(p) {
  a <- gene_models[gene_group == p[1], adj_r2]
  b2 <- gene_models[gene_group == p[2], adj_r2]
  obs <- median(a) - median(b2)
  bs  <- replicate(n_boot,
                   median(sample(a, replace = TRUE)) - median(sample(b2, replace = TRUE)))
  ci  <- unname(quantile(bs, c(0.025, 0.975)))
  data.table(class_a = p[1], class_b = p[2], n_a = length(a), n_b = length(b2),
             diff_median_adj_r2 = obs, ci_lo = ci[1], ci_hi = ci[2],
             equivalent = abs(ci[1]) < r2_margin & abs(ci[2]) < r2_margin)
}))
fwrite(equiv, outf("class_equivalence", "csv"))

cat(sprintf("\n== 8b. Between-class equivalence (margin +/- %.2f adjusted R2) ==\n", r2_margin))
print(equiv[, .(class_a, class_b, n_a, n_b,
                diff = round(diff_median_adj_r2, 3),
                CI = paste0(round(ci_lo, 3), " to ", round(ci_hi, 3)),
                equivalent)], row.names = FALSE)
if (all(equiv$equivalent)) {
  cat(sprintf("All pairwise CIs fall within +/- %.2f: equivalent at this margin, which\n", r2_margin))
  cat("is a stronger statement than non-significance.\n")
} else {
  cat(sprintf("At least one CI extends beyond +/- %.2f, so equivalence is NOT\n", r2_margin))
  cat("established at this margin. Report the interval, not the p-value.\n")
}
cat("CAVEAT: the bootstrap resamples genes independently. Co-expressed genes are\n")
cat("not independent, so these intervals are narrower than the truth and this\n")
cat("test errs toward declaring equivalence. Read it as an upper bound on\n")
cat("precision, and do not rest a conclusion on a CI that only just fits.\n")

## ============================================================================
## 8c -- do the classes differ in WHICH cell type they track?
## ============================================================================
## The artefact hypothesis is directional: if composition drives the labels,
## reference-enriched genes should be neuron-associated and case-enriched genes
## glia-associated. An omnibus test cannot see direction, so the class pair is
## compared one-sided where the two enriched classes can be identified.
ce <- if (!is.null(params$case_enriched_group)) params$case_enriched_group else
        grep(case, classes, value = TRUE, fixed = TRUE)[1]
re <- if (!is.null(params$ref_enriched_group)) params$ref_enriched_group else
        grep(ctrl, classes, value = TRUE, fixed = TRUE)[1]
directional <- !is.na(ce) && !is.na(re) && ce %in% classes && re %in% classes

class_cmp <- rbindlist(lapply(keep_ct, function(ct) {
  col <- paste0("beta_", ct)
  k <- suppressWarnings(kruskal.test(gene_models[[col]], factor(gene_models$gene_group)))
  med <- gene_models[, .(m = median(get(col))), by = gene_group]
  row <- data.table(cell_type = ct,
                    medians = paste(sprintf("%s=%.3f", med$gene_group, med$m), collapse = "; "),
                    pval_omnibus = k$p.value,
                    pval_directional = NA_real_, direction = NA_character_)
  if (directional) {
    x <- gene_models[gene_group == ce, get(col)]
    y <- gene_models[gene_group == re, get(col)]
    if (length(x) >= 3L && length(y) >= 3L) {
      ## Two one-sided tests, both reported, so the direction of any difference
      ## is explicit rather than inferred from a two-sided p-value.
      hi <- suppressWarnings(wilcox.test(x, y, alternative = "greater", exact = FALSE))
      lo <- suppressWarnings(wilcox.test(x, y, alternative = "less",    exact = FALSE))
      row[, `:=`(pval_directional = min(hi$p.value, lo$p.value),
                 direction = if (hi$p.value < lo$p.value)
                               sprintf("%s > %s", ce, re) else sprintf("%s < %s", ce, re))]
    }
  }
  row
}))
class_cmp[, fdr_omnibus := p.adjust(pval_omnibus, "BH")]
if (directional) class_cmp[, fdr_directional := p.adjust(pval_directional, "BH")]
fwrite(class_cmp, outf("class_comparison", "csv"))

cat("\n== 8c. Cell-type association by class ==\n")
if (!directional) {
  cat("Enriched-class labels could not be identified from the group names, so only\n")
  cat("the omnibus test is reported. Set --case-enriched-group/--ref-enriched-group\n")
  cat("to get the directional test the artefact hypothesis actually predicts.\n")
}
print(class_cmp, row.names = FALSE)

## ============================================================================
## 8f -- THE DIRECT TEST: does composition account for the CN/AD call difference?
## ============================================================================
## Everything above is indirect. This re-estimates the group difference in the
## CALLS themselves with composition in the model, for every labelled gene.
##
## Outcome: whether an individual carries any call in the gene -- the same
## presence/absence the global labels are built from. Three nested models, so
## the movement in the diagnosis coefficient is visible rather than asserted:
##
##   M1  present ~ diagnosis                       the uncorrected contrast
##   M2  M1 + cell-type proportions                composition adjusted
##   M3  M2 + log2 expression + log2 library size  detection adjusted
##
## M3 is the estimate to report. Expression and library size are the observed
## determinants of depth, and adjusting for them conditions on what makes a
## position testable WITHOUT imputing which zeros were untestable. Because each
## model is fitted within one gene, sample-constant coverage biases are absorbed.
##
## Complete separation is not an edge case: a gene called only in cases is
## exactly the presence/absence finding the mechanism is alleged to manufacture,
## and it drives the ML estimate to infinity. Firth's penalised likelihood is the
## right estimator and is used when logistf is installed.
has_firth <- requireNamespace("logistf", quietly = TRUE)
if (!has_firth) {
  cat("\nNOTE: package 'logistf' is not installed, so genes whose calls separate\n")
  cat("perfectly on diagnosis are flagged and excluded from 8f rather than fitted.\n")
  cat("Those are the presence/absence genes -- the most relevant cases. Install\n")
  cat("logistf to include them.\n")
}

carriers <- unique(sites[, .(ref_base, rids)])[, present := 1L]
call_tab <- CJ(ref_base = unique(expr_long$ref_base), rids = wide$rids, unique = TRUE)
call_tab <- merge(call_tab, carriers, by = c("ref_base", "rids"), all.x = TRUE)
call_tab[is.na(present), present := 0L]
call_tab <- merge(call_tab, expr_long, by = c("ref_base", "rids"))
call_tab <- merge(call_tab, wide, by = "rids")
call_tab[, `:=`(diagnosis = factor(con, levels = c(ctrl, case)),
                lexpr = log2(expression + 1),
                llib  = log2(lib_size[rids]))]

gname  <- paste0("diagnosis", case)
## Braced deliberately: Rscript evaluates top-level expressions one at a time,
## so an unbraced `if` whose `else` starts the next line parses as a complete
## expression followed by a stray `else`.
fitter <- if (has_firth) {
  function(form, dat) tryCatch(logistf::logistf(form, data = dat), error = function(e) NULL)
} else {
  function(form, dat) tryCatch(glm(form, family = binomial, data = dat), error = function(e) NULL)
}

grab <- function(f) {
  if (is.null(f)) return(c(NA_real_, NA_real_))
  if (inherits(f, "logistf")) {
    i <- match(gname, names(coef(f)))
    return(if (is.na(i)) c(NA_real_, NA_real_) else c(unname(coef(f)[i]), unname(f$prob[i])))
  }
  cf <- coef(summary(f))
  if (!gname %in% rownames(cf)) return(c(NA_real_, NA_real_))
  c(cf[gname, "Estimate"], cf[gname, "Pr(>|z|)"])
}
## Firth is built for separation, so only the ML path needs the check.
is_sep <- function(f) {
  if (is.null(f) || inherits(f, "logistf")) return(FALSE)
  any(fitted(f) < 1e-6 | fitted(f) > 1 - 1e-6) || any(abs(coef(f)) > 20, na.rm = TRUE)
}

f_comp <- as.formula(paste("present ~ diagnosis +", paste(keep_ct, collapse = " + ")))
f_det  <- as.formula(paste("present ~ diagnosis +", paste(keep_ct, collapse = " + "),
                           "+ lexpr + llib"))

fit_one_gene <- function(s) {
  n_pos <- sum(s$present); n_neg <- nrow(s) - n_pos
  if (n_pos < min_carriers || n_neg < min_carriers) return(NULL)
  m1 <- fitter(present ~ diagnosis, s)
  m2 <- fitter(f_comp, s)
  m3 <- if (sd(s$lexpr) > 0) fitter(f_det, s) else NULL
  a <- grab(m1); b2 <- grab(m2); c2 <- grab(m3)
  data.table(ref_base = s$ref_base[1], n = nrow(s), n_carriers = n_pos,
             n_carriers_case = sum(s$present == 1L & s$con == case),
             n_carriers_ref  = sum(s$present == 1L & s$con == ctrl),
             min_expression  = min(s$expression),
             beta_naive = a[1],  p_naive = a[2],
             beta_comp  = b2[1], p_comp  = b2[2],
             beta_det   = c2[1], p_det   = c2[2],
             separated  = is_sep(m1) || is_sep(m2) || is_sep(m3),
             estimator  = if (has_firth) "firth" else "ml")
}

call_models <- rbindlist(lapply(split(call_tab, by = "ref_base"), fit_one_gene), fill = TRUE)

if (nrow(call_models)) {
  att <- function(b0, b1) fifelse(is.na(b0) | abs(b0) < 1e-3, NA_real_,
                                  100 * (1 - abs(b1) / abs(b0)))
  call_models[, `:=`(att_comp = att(beta_naive, beta_comp),
                     att_det  = att(beta_naive, beta_det))]
  call_models[separated == FALSE, `:=`(fdr_naive = p.adjust(p_naive, "BH"),
                                       fdr_comp  = p.adjust(p_comp,  "BH"),
                                       fdr_det   = p.adjust(p_det,   "BH"))]
  call_models <- merge(call_models, labels, by = "ref_base", all.x = TRUE)
  fwrite(call_models, outf("call_models", "csv"))

  ok  <- call_models[separated == FALSE]
  sig <- ok[!is.na(fdr_naive) & fdr_naive < 0.05]

  cat("\n== 8f. Gene-level call presence: naive vs composition vs detection ==\n")
  cat(sprintf("%d gene(s) fitted (>= %d carriers and non-carriers); %d excluded for\n",
              nrow(call_models), min_carriers, sum(call_models$separated)))
  cat(sprintf("separation. Estimator: %s.\n", paste(unique(call_models$estimator), collapse = "/")))
  cat(sprintf("Significant at FDR < 0.05 -- naive %d, composition-adjusted %d, \n",
              nrow(sig), sum(ok$fdr_comp < 0.05, na.rm = TRUE)))
  cat(sprintf("detection-adjusted %d.\n", sum(ok$fdr_det < 0.05, na.rm = TRUE)))

  ## Reported over ALL fitted genes as well as the significant ones. Conditioning
  ## an attenuation summary on a significant naive effect selects on the very
  ## quantity in its denominator and inflates the figure.
  cat(sprintf("\nMedian attenuation over all %d fitted genes: composition %.1f%%, detection %.1f%%\n",
              nrow(ok), median(ok$att_comp, na.rm = TRUE), median(ok$att_det, na.rm = TRUE)))
  cat(sprintf("Median attenuation among the %d significant naively: composition %.1f%%, detection %.1f%%\n",
              nrow(sig), median(sig$att_comp, na.rm = TRUE), median(sig$att_det, na.rm = TRUE)))
  cat("The first line is the honest estimate; the second is selected on the naive\n")
  cat("effect and is given only because it is the number a reader looks for.\n")

  cat("\nBy class:\n")
  print(ok[, .(n = .N,
               med_beta_naive = round(median(beta_naive, na.rm = TRUE), 3),
               med_beta_comp  = round(median(beta_comp,  na.rm = TRUE), 3),
               med_beta_det   = round(median(beta_det,   na.rm = TRUE), 3),
               med_att_det    = round(median(att_det,    na.rm = TRUE), 1),
               n_sig_naive = sum(fdr_naive < 0.05, na.rm = TRUE),
               n_sig_det   = sum(fdr_det   < 0.05, na.rm = TRUE)),
           by = gene_group], row.names = FALSE)

  ## ---- sensitivity to the risk set ------------------------------------------
  ## A zero is ambiguous only where the gene might not have been testable. For
  ## genes expressed well above the floor in EVERY sample the ambiguity does not
  ## arise, so the analysis is repeated over a range of minimum-expression cut-
  ## offs. No single cut-off is asserted: the point is whether the answer is
  ## stable as the risk set is tightened. If it is, the ambiguous zeros are not
  ## driving it. If it moves, that is the finding.
  cutoffs <- unique(round(quantile(call_models$min_expression,
                                   c(0, .25, .5, .75), na.rm = TRUE), 4))
  sens <- rbindlist(lapply(cutoffs, function(cut) {
    s <- ok[min_expression >= cut]
    if (!nrow(s)) return(NULL)
    data.table(min_expression_cutoff = cut, n_genes = nrow(s),
               median_att_comp = median(s$att_comp, na.rm = TRUE),
               median_att_det  = median(s$att_det,  na.rm = TRUE),
               n_sig_naive = sum(s$fdr_naive < 0.05, na.rm = TRUE),
               n_sig_det   = sum(s$fdr_det   < 0.05, na.rm = TRUE))
  }))
  fwrite(sens, outf("call_sensitivity", "csv"))
  cat("\n== 8f. Sensitivity to the risk set (minimum expression across samples) ==\n")
  print(sens[, .(cutoff = min_expression_cutoff, n_genes,
                 att_comp = round(median_att_comp, 1),
                 att_det  = round(median_att_det, 1),
                 n_sig_naive, n_sig_det)], row.names = FALSE)
  cat("Rows to the right of the first tighten the risk set to genes whose zeros are\n")
  cat("progressively less ambiguous. A stable column means the ambiguity is not\n")
  cat("carrying the result.\n")
}

## ============================================================================
## 8g -- the middle link, as far as it is observable
## ============================================================================
## The mechanism is composition -> expression -> depth -> threshold -> call. 8b
## measures the first arrow. This measures the second, using the only depth this
## data has.
##
## LIMITATION, and it is not small: dp is recorded only where a call was made,
## so this is depth CONDITIONAL ON DETECTION. It says nothing about depth at
## positions that never cleared the threshold -- the positions the mechanism is
## about. It is a lower bound on the composition-depth association and it feeds
## nothing downstream; no value here is extrapolated anywhere.
lvl_file <- file.path(out_dir, "08_site_level.csv")
lvl <- if (file.exists(lvl_file)) fread(lvl_file)[dp > 0] else NULL

if (!is.null(lvl) && nrow(lvl)) {
  dpg <- lvl[, .(median_dp = median(dp)), by = .(ref_base, rids)]
  dpg <- merge(dpg, wide, by = "rids")
  dpg[, ldp := log2(median_dp + 1)]
  dp_form <- as.formula(paste("ldp ~", paste(keep_ct, collapse = " + ")))
  depth_models <- rbindlist(lapply(split(dpg, by = "ref_base"), function(s) {
    if (nrow(s) < min_n || sd(s$ldp) == 0) return(NULL)
    fit <- lm(dp_form, data = s)
    data.table(ref_base = s$ref_base[1], n = nrow(s),
               r2 = summary(fit)$r.squared, adj_r2 = summary(fit)$adj.r.squared,
               pval = tryCatch({ f <- summary(fit)$fstatistic
                                 unname(pf(f[1], f[2], f[3], lower.tail = FALSE)) },
                               error = function(e) NA_real_))
  }), fill = TRUE)
  if (nrow(depth_models)) {
    depth_models[, fdr := p.adjust(pval, "BH")]
    depth_models <- merge(depth_models, labels, by = "ref_base", all.x = TRUE)
    fwrite(depth_models, outf("depth_models", "csv"))
    cat("\n== 8g. Composition R2 for observed read depth, by class ==\n")
    print(depth_models[, .(n = .N, median_adj_r2 = round(median(adj_r2), 3),
                           pct_fdr05 = round(100 * mean(fdr < 0.05), 1)),
                       by = gene_group], row.names = FALSE)
    cat("Depth here is conditional on a call having been made. A high value means\n")
    cat("the composition-to-depth arrow is real for that gene; a low one does NOT\n")
    cat("clear it, because uncalled positions are invisible to this model.\n")
  }
} else {
  depth_models <- NULL
}

## ============================================================================
## 8d -- site-level adjustment, every labelled gene
## ============================================================================
## No longer gated on --site-genes: this and 8f are the only places a CN/AD
## effect is re-estimated with composition in it, and gating that on a
## hand-typed list left the central question answered for a handful of genes.
##
## The outcome is a binomial proportion whose variance depends on dp, so OLS
## would weight a position at dp=20 like one at dp=400. Quasibinomial keeps the
## mean model and lets the data set the dispersion.
##
## The old depth-adjusted model regressed a1ad/dp on log(dp) -- the outcome's own
## denominator, which manufactures an association whatever the biology. Depth is
## tested where it belongs instead, as an outcome, in 8g.
site_models <- NULL
if (!is.null(lvl) && nrow(lvl)) {
  lvl2 <- merge(lvl, wide[, c("rids", keep_ct), with = FALSE], by = "rids")
  lvl2[, diagnosis := factor(con, levels = c(ctrl, case))]
  s_form <- as.formula(paste("cbind(a1ad, dp - a1ad) ~ diagnosis +",
                             paste(keep_ct, collapse = " + ")))

  ## Positions are filtered BEFORE the split, not inside the loop. Now that
  ## every labelled gene is loaded this is the difference between fitting two
  ## GLMs per position for tens of thousands of positions and fitting them only
  ## where the result is interpretable. A model with one parameter per retained
  ## cell type plus intercept and diagnosis needs the residual degrees of
  ## freedom, so --min-n is the cut, and positions below it are not fitted at
  ## all rather than fitted and then discarded.
  lvl2[, `:=`(n_pos = .N, n_grp = uniqueN(diagnosis)), by = .(ref_base, poss, a1type)]
  n_positions_seen <- uniqueN(lvl2[, .(ref_base, poss, a1type)])
  lvl2 <- lvl2[n_pos >= min_n & n_grp >= 2L]
  n_positions_fit <- uniqueN(lvl2[, .(ref_base, poss, a1type)])
  cat(sprintf("\n8d: %d of %d position(s) have n >= %d in both groups and are fitted.\n",
              n_positions_fit, n_positions_seen, min_n))

  site_models <- if (!nrow(lvl2)) NULL else
    rbindlist(lapply(split(lvl2, by = c("ref_base", "poss", "a1type"), drop = TRUE), function(s) {
    ## Warnings are suppressed because the condition that raises them -- fitted
    ## probabilities at 0 or 1 -- is detected explicitly below and acted on.
    ## Left unsuppressed it is one warning per boundary position, which buries
    ## everything else the step prints.
    f0 <- suppressWarnings(tryCatch(glm(cbind(a1ad, dp - a1ad) ~ diagnosis,
                       family = quasibinomial, data = s), error = function(e) NULL))
    f1 <- suppressWarnings(tryCatch(glm(s_form, family = quasibinomial, data = s),
                   error = function(e) NULL))
    g <- function(f) {
      if (is.null(f)) return(c(NA_real_, NA_real_))
      cf <- coef(summary(f)); nm <- paste0("diagnosis", case)
      if (!nm %in% rownames(cf)) return(c(NA_real_, NA_real_))
      c(cf[nm, "Estimate"], cf[nm, "Pr(>|t|)"])
    }
    ## A position where one group is entirely unedited (or entirely edited)
    ## drives the diagnosis coefficient to the boundary: glm returns a huge
    ## estimate with a meaningless standard error. That is not a test, so it is
    ## flagged and kept out of the BH family -- the same treatment 8f gives
    ## separated genes. The rows stay in the CSV with the flag.
    at_bound <- function(f) {
      if (is.null(f)) return(TRUE)
      !isTRUE(f$converged) ||
        any(fitted(f) < 1e-8 | fitted(f) > 1 - 1e-8) ||
        any(abs(coef(f)) > 10, na.rm = TRUE)
    }
    a <- g(f0); b2 <- g(f1)
    data.table(ref_base = s$ref_base[1], poss = s$poss[1], a1type = s$a1type[1],
               n = nrow(s), boundary = at_bound(f0) || at_bound(f1),
               beta_naive = a[1], p_naive = a[2],
               beta_comp = b2[1], p_comp = b2[2],
               pct_attenuation = if (is.na(a[1]) || abs(a[1]) < 1e-3) NA_real_
                                 else 100 * (1 - abs(b2[1]) / abs(a[1])))
  }), fill = TRUE)

  if (!is.null(site_models) && nrow(site_models)) {
    ## Every position here cleared --min-n, so the only exclusion left is the
    ## boundary one: an unbounded coefficient is not a test, and leaving it in
    ## the family would deflate the correction for every position that is.
    site_models[boundary == FALSE, `:=`(fdr_naive = p.adjust(p_naive, "BH"),
                                        fdr_comp  = p.adjust(p_comp,  "BH"))]
    fwrite(site_models, outf("site_models", "csv"))

    est <- site_models[boundary == FALSE]
    sig <- est[!is.na(fdr_naive) & fdr_naive < 0.05]
    cat("\n== 8d. Site-level adjustment (quasibinomial, all labelled genes) ==\n")
    cat(sprintf("%d position(s) modelled at n >= %d; %d at the boundary (one group\n",
                nrow(site_models), min_n, sum(site_models$boundary)))
    cat(sprintf("uniformly unedited or edited) and excluded, leaving %d.\n", nrow(est)))
    cat(sprintf("%d significant before adjustment (FDR < 0.05); %d after adjusting\n",
                nrow(sig), sum(sig$fdr_comp < 0.05, na.rm = TRUE)))
    cat("for composition.\n")
    cat(sprintf("Median attenuation over all %d positions: %.1f%%\n",
                nrow(est), median(est$pct_attenuation, na.rm = TRUE)))
    if (nrow(sig)) {
      cat(sprintf("Median attenuation among the %d significant naively: %.1f%%\n",
                  nrow(sig), median(sig$pct_attenuation, na.rm = TRUE)))
      cat("\nTop positions significant before adjustment:\n")
      print(head(sig[order(fdr_naive), .(ref_base, poss, a1type, n,
                                         beta_naive = signif(beta_naive, 3),
                                         beta_comp = signif(beta_comp, 3),
                                         pct_att = round(pct_attenuation, 1),
                                         fdr_naive = signif(fdr_naive, 3),
                                         fdr_comp = signif(fdr_comp, 3))], 25L),
            row.names = FALSE)
    }
    cat("\nPositions called in only one group cannot support this model -- there is\n")
    cat("nothing to compare. Whether the other group was even testable at those\n")
    cat("positions is NOT recoverable from this database; it needs the VCFs or BAMs.\n")
    cat("No estimate of it is offered here.\n")
  }
}

## ============================================================================
## 8e -- named genes located within the full distribution
## ============================================================================
site_genes <- params$site_genes
if (length(site_genes)) {
  sub <- gene_models[gene_name %in% site_genes]
  if (nrow(sub)) {
    sub[, r2_percentile := vapply(adj_r2, function(x) mean(gene_models$adj_r2 <= x), numeric(1))]
    cat("\n== 8e. Named genes against the full distribution ==\n")
    print(sub[order(-adj_r2), .(gene_name, gene_group, adj_r2 = round(adj_r2, 3),
                                percentile = round(100 * r2_percentile, 1),
                                fdr = signif(fdr, 3))], row.names = FALSE)
  }
  keys <- names(sym_of)[sym_of %in% site_genes]
  if (nrow(call_models) && "separated" %in% names(call_models)) {
    cm <- call_models[ref_base %in% keys]
    if (nrow(cm)) {
      cat("\n8f result for the named genes:\n")
      print(cm[, .(ref_base, gene_group, n_carriers_case, n_carriers_ref,
                   beta_naive = signif(beta_naive, 3),
                   beta_det = signif(beta_det, 3),
                   att_det = round(att_det, 1), separated)], row.names = FALSE)
    }
  }
}

## ---- what this analysis cannot settle --------------------------------------
cat("\n== Limits ==\n")
cat("A (gene, sample) pair with no call is either a true negative or a position\n")
cat("that never had the depth to be tested. Which of the two it is cannot be\n")
cat("recovered from this database: dp exists only where a call was made, and the\n")
cat("VCFs and BAMs are gone. Every model above adjusts for the OBSERVED drivers\n")
cat("of depth -- gene expression and library size -- rather than guessing which\n")
cat("zeros were untestable, and each is fitted within a gene so sample-constant\n")
cat("coverage bias is absorbed. What remains uncontrolled is sample-by-gene\n")
cat("interaction in coverage bias. There is no basis in this data for a single\n")
cat("attributable-fraction figure, and none is reported.\n")

## ---- plot ------------------------------------------------------------------
png(outf("confounding", "png"), width = 1800, height = 700, res = 150)
par(mfrow = c(1, 2))
boxplot(adj_r2 ~ gene_group, data = gene_models, ylab = "composition adj R2",
        xlab = "", main = "8b: composition R2 by class", col = "steelblue")
if (nrow(call_models) && "separated" %in% names(call_models)) {
  ## Composition-adjusted, not detection-adjusted: this panel supports the
  ## reported result (composition changes nothing). The detection-adjusted
  ## coefficient is in 08_call_models.csv for anyone who wants it.
  okp <- call_models[separated == FALSE & !is.na(beta_naive) & !is.na(beta_comp)]
  plot(okp$beta_naive, okp$beta_comp, pch = 19, cex = 0.5,
       col = adjustcolor("black", 0.5),
       xlab = "diagnosis log-odds, naive", ylab = "diagnosis log-odds, composition-adjusted",
       main = "8f: movement under adjustment")
  abline(0, 1, col = "grey50"); abline(h = 0, v = 0, lty = 3)
}
invisible(dev.off())

write_json(list(label = params$label, case_group = case, reference_group = ctrl,
                aggregator = agg, expression = all_name,
                dropped_cell_type = drop_ct, seed = seed,
                dp_threshold = dp_thresh, min_carriers = min_carriers,
                burden_ratio = unname(per_capita[["case"]] / per_capita[["ref"]]),
                burden_pval = bw$p.value, burden_vs_composition = burden_cmp,
                n_genes_modelled = nrow(gene_models),
                class_r2 = gene_models[, .(n = .N, median_adj_r2 = median(adj_r2)),
                                       by = gene_group],
                kruskal_adj_r2_p = kw$p.value, class_comparison = class_cmp,
                r2_equivalence = list(margin = r2_margin, n_bootstrap = n_boot,
                                      pairs = equiv,
                                      all_equivalent = all(equiv$equivalent),
                                      assumes_gene_independence = TRUE),
                call_models = if (nrow(call_models)) call_models else NULL,
                call_sensitivity = if (exists("sens")) sens else NULL,
                depth_models_conditional_on_detection = depth_models,
                site_min_n = min_n),
           outf("confounding_summary", "json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("\nWrote results to %s\n", out_dir))
