#!/usr/bin/env Rscript
##
## 02_build_matrix.R -- gene-level expression from StringTie ctab output.
##
##   Rscript 02_build_matrix.R <params.json>
##
## Transcripts reach a symbol through gene_id and the step-1 locus map, never
## through the ctab gene_name column, which is "." for novel isoforms inside
## annotated loci -- 35.7% of transcripts in this dataset.
##

suppressPackageStartupMessages({
  library(jsonlite); library(data.table); library(HGNChelper)
  library(edgeR); library(BRETIGEA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 02_build_matrix.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir   <- params$output_dir
read_len  <- params$read_length
min_fpkm  <- params$min_marker_fpkm
enst_only <- isTRUE(params$enst_only)

## Every row is a matched cohort sample -- Python refuses to write this file
## otherwise -- so no filtering is needed here.
crosswalk <- fread(file.path(out_dir, "02_sample_crosswalk.csv"))
if (nrow(crosswalk) == 0L) stop("no samples in 02_sample_crosswalk.csv")
## Step 1's map contains clean loci only, so the join below drops chimeric loci
## by itself. The locus list is read purely to attribute the loss in the report.
gmap      <- fread(file.path(out_dir, "01_locus_symbol_map.csv"))
chim_loci <- readLines(file.path(out_dir, "01_chimeric_loci.txt"))

## ---- ingest ----------------------------------------------------------------
td <- rbindlist(Map(function(f, s) {
  fread(f, select = c("t_name","gene_id","gene_name","length","cov","FPKM"),
        showProgress = FALSE)[, sample := s]
}, crosswalk$ctab_path, crosswalk$rids))
n_tx_total <- nrow(td)

## ---- symbols, two routes in priority order ---------------------------------
## 1. the ctab's own gene_name -- the annotation stating outright which gene a
##    transcript belongs to. Unambiguous, so chimerism is irrelevant here:
##    SNAP25's transcripts are SNAP25's whether or not StringTie fused their
##    locus with a neighbour.
## 2. for transcripts with no gene_name (novel isoforms), the locus map --
##    inherit the symbol when the locus resolves to exactly one gene, drop when
##    it is chimeric or carries no reference gene at all.
##
## Resolving through the locus map alone, as an earlier version did, discarded
## named transcripts sitting in fused loci. That removed SNAP25, RBFOX3, MBP and
## CLDN5 -- the highest-expressed neuronal, oligodendrocyte and endothelial
## markers -- while leaving GFAP and AQP4 intact, which handed astrocytes half
## the estimated mRNA pool.
td[, gene_name := fifelse(gene_name %in% c(".", ""), NA_character_, gene_name)]
td <- merge(td, gmap[!duplicated(locus)], by.x = "gene_id", by.y = "locus", all.x = TRUE)
td[, symbol := fifelse(!is.na(gene_name), gene_name, symbol)]

n_named     <- sum(!is.na(td$gene_name))
n_inherited <- sum(is.na(td$gene_name) & !is.na(td$symbol))
n_chim_rows <- sum(is.na(td$symbol) & td$gene_id %in% chim_loci)
n_unnamed   <- sum(is.na(td$symbol)) - n_chim_rows
td <- td[!is.na(symbol)]
if (enst_only) td <- td[grepl("^ENST", t_name)]
if (nrow(td) == 0L) stop("no transcripts survived symbol resolution and filtering")

## ---- gene-level matrices ---------------------------------------------------
gfpkm <- dcast(td, symbol ~ sample, value.var = "FPKM", fun.aggregate = sum, fill = 0)
mat_fpkm <- as.matrix(gfpkm[, -1]); rownames(mat_fpkm) <- gfpkm$symbol

td[, est_counts := cov * length / read_len]
gcnt <- dcast(td, symbol ~ sample, value.var = "est_counts", fun.aggregate = sum, fill = 0)
mat_cnt <- round(as.matrix(gcnt[, -1])); rownames(mat_cnt) <- gcnt$symbol
mat_cnt <- mat_cnt[rownames(mat_fpkm), colnames(mat_fpkm), drop = FALSE]

## ---- harmonise symbols to current HGNC -------------------------------------
## BRETIGEA's lists are 2018-vintage; without this, aliased symbols silently
## shrink the marker sets.
fix <- checkGeneSymbols(rownames(mat_fpkm), species = "human")
new <- ifelse(is.na(fix$Suggested.Symbol), rownames(mat_fpkm), fix$Suggested.Symbol)
new <- sub(" ///.*$", "", new)
n_sym_before <- nrow(mat_fpkm)
old_names <- rownames(mat_fpkm)
mat_fpkm <- rowsum(mat_fpkm, group = new)
mat_cnt  <- rowsum(mat_cnt,  group = new)

write_mat <- function(m, f) fwrite(data.table(symbol = rownames(m), as.data.table(m)),
                                   file.path(out_dir, f))

## ---- exclude by annotation before normalising ------------------------------
## The gene list comes from step 1's --biotype-tsv: mitochondrial, rRNA, tRNA
## and misc_RNA by default, identified by what the annotation says they are
## rather than by their symbols. Their abundance reflects rRNA-depletion
## efficiency and RNA quality, not cell content -- and because TPM forces a
## fixed per-sample total, their variation compresses every real gene by a
## different amount in every sample. That sample-varying compression becomes
## the dominant axis of the matrix and swamps the composition signal
## deconvolution depends on. Removing them here means the denominator is
## actual nuclear mRNA.
excl_ids <- character(0)
gene_file <- file.path(out_dir, "00_excluded_genes.csv")
if (file.exists(gene_file)) {
  ex <- fread(gene_file, colClasses = "character")
  excl_ids <- c(ex$gene_id, ex$gene_id_base)
}
sym_file <- file.path(out_dir, "01_excluded_symbols.txt")
if (file.exists(sym_file)) excl_ids <- c(excl_ids, readLines(sym_file))
excl_ids <- unique(excl_ids[!is.na(excl_ids) & nzchar(excl_ids)])

## A harmonised symbol is excluded when it is on the list itself or when any
## symbol folded into it was -- HGNChelper can rename a gene out of the list.
excl_new <- unique(c(intersect(rownames(mat_fpkm), excl_ids), new[old_names %in% excl_ids]))
excl <- rownames(mat_fpkm) %in% excl_new

lib        <- colSums(mat_fpkm)
excl_frac  <- colSums(mat_fpkm[excl, , drop = FALSE]) / lib
n_excluded <- sum(excl)

## Read by step 3, which correlates it against the leading components. If PC1
## still tracks the fraction removed, the leading axis of the matrix is still
## library composition rather than biology.
fwrite(data.table(rids = colnames(mat_fpkm), fpkm_total = lib,
                  frac_excluded = excl_frac),
       file.path(out_dir, "02_excluded_fraction.csv"))

mat_fpkm <- mat_fpkm[!excl, , drop = FALSE]
mat_cnt  <- mat_cnt[!excl, , drop = FALSE]

## ---- what still dominates --------------------------------------------------
## The check on whether the exclusion was sufficient. Any single gene holding a
## percent or more of the retained budget is doing to the TPM denominator what
## the structural ncRNA did, and is worth looking up in the annotation table.
top_share <- head(sort(rowMeans(sweep(mat_fpkm, 2, colSums(mat_fpkm), "/")),
                       decreasing = TRUE), 12L)

## ---- filter, rescuing markers ----------------------------------------------
grp <- crosswalk$con[match(colnames(mat_cnt), crosswalk$rids)]
keep <- filterByExpr(DGEList(mat_cnt), group = grp)
marker_universe <- unique(BRETIGEA::markers_df_human_brain$markers)

## The rescue threshold must be applied on a normalised scale. StringTie's FPKM
## was computed against the full library including the structural ncRNA
## fraction, which varied 50-92% between samples -- so a fixed cut on raw FPKM
## means something different in each sample, and is systematically stricter in
## the samples where that fraction was largest. Rescaling to the retained
## denominator first makes the threshold comparable across samples.
fpkm_norm <- sweep(mat_fpkm, 2, colSums(mat_fpkm), "/") * 1e6
rescued <- rownames(mat_fpkm) %in% marker_universe & rowMeans(fpkm_norm) > min_fpkm
n_rescued_raw  <- sum(rownames(mat_fpkm) %in% marker_universe &
                      rowMeans(mat_fpkm) > min_fpkm & !keep)
n_rescued <- sum(rescued & !keep)
keep <- keep | rescued

## Unfiltered matrices are set aside before the subset. Step 8 models labelled
## genes filterByExpr removed -- disproportionately the low-expressed ones
## nearest the dp >= 20 boundary, which are exactly the genes most exposed to
## the detection mechanism it is testing.
fpkm_all    <- mat_fpkm
cnt_all     <- mat_cnt
symbols_all <- rownames(mat_fpkm)

mat_fpkm <- mat_fpkm[keep, , drop = FALSE]
mat_cnt  <- mat_cnt[keep, , drop = FALSE]

## ---- TPM, linear scale for DSA ---------------------------------------------
## Computed after filtering, so the denominator is the retained set. Per-sample
## rescaling is exactly what TPM does; DSA needs linear input and renormalises
## proportions itself, so this is the intended scale.
tpm <- sweep(mat_fpkm, 2, colSums(mat_fpkm), "/") * 1e6

## ---- TMM as an alternative linear normalisation -----------------------------
## TMM derives each sample's scaling from a trimmed mean of log-ratios, so a
## minority of dominant genes cannot move it -- including high-expressing loci
## with no symbol, which no exclusion pattern can catch. cpm() returns linear
## values, so this remains valid input for DSA. Step 5 chooses between the two
## with --expression.
dge <- calcNormFactors(DGEList(mat_cnt), method = "TMM")
tmm <- cpm(dge, normalized.lib.sizes = TRUE)

## ---- unfiltered matrices, on the FILTERED sample scale ----------------------
## Step 8 needs expression for whatever gene set the editing analysis names,
## which filterByExpr has no reason to have retained. The denominators here are
## the filtered ones, applied to the unfiltered rows. A different per-sample
## denominator is a per-sample additive shift on log2, and since the estimated
## proportions are themselves derived per sample from the filtered matrix, a
## second scale can induce or mask exactly the association step 8 exists to
## measure. On the filtered scale these are a strict superset: identical values
## on shared genes, extra rows for the rest.
##
## TMM factors are estimated on the filtered matrix on purpose -- edgeR's
## trimmed mean is unreliable on a matrix that is mostly near-zero -- then
## carried over rather than re-estimated, so the scaling is the same one the
## deconvolution ran on.
tpm_all <- sweep(fpkm_all, 2, colSums(mat_fpkm), "/") * 1e6
tmm_all <- cpm(DGEList(counts       = cnt_all,
                       lib.size     = dge$samples$lib.size,
                       norm.factors = dge$samples$norm.factors),
               normalized.lib.sizes = TRUE)

## The superset property is the whole point of the construction, so it is
## asserted rather than assumed -- a future change to either denominator would
## otherwise reintroduce the two-scale problem silently.
stopifnot(all.equal(tpm_all[rownames(tpm), , drop = FALSE], tpm),
          all.equal(tmm_all[rownames(tmm), , drop = FALSE], tmm))

## ---- marker survival -------------------------------------------------------
## The number that actually governs whether step 5 is well-posed: step 1's
## counts were measured in the GTF, before HGNC harmonisation and before any
## expression filtering, so they do not describe this matrix.
md    <- BRETIGEA::markers_df_human_brain
top_n <- if (is.null(params$top_n)) 500L else params$top_n
marker_survival <- rbindlist(lapply(sort(unique(md$cell)), function(ct) {
  top <- head(md$markers[md$cell == ct], top_n)
  data.table(cell_type = ct,
             in_matrix_all      = length(intersect(top, symbols_all)),
             in_matrix_filtered = length(intersect(top, rownames(tpm))))
}))
fwrite(marker_survival, file.path(out_dir, "02_marker_survival.csv"))

## ---- outputs ---------------------------------------------------------------
write_mat(mat_fpkm, "02_gene_fpkm.csv")
write_mat(mat_cnt,  "02_gene_counts.csv")
write_mat(tpm,      "02_gene_tpm.csv")
write_mat(tmm,      "02_gene_tmm.csv")
write_mat(tpm_all,  "02_gene_tpm_all.csv")
write_mat(tmm_all,  "02_gene_tmm_all.csv")
saveRDS(list(fpkm = mat_fpkm, counts = mat_cnt, tpm = tpm, tmm = tmm,
             tpm_all = tpm_all, tmm_all = tmm_all, group = grp),
        file.path(out_dir, "02_expression.rds"))

summary_list <- list(
  n_samples            = ncol(tpm),
  read_length          = read_len,
  enst_only            = enst_only,
  n_transcripts_read   = n_tx_total,
  n_tx_named_directly  = n_named,
  n_tx_symbol_inherited = n_inherited,
  n_tx_unresolved      = n_unnamed,
  n_tx_chimeric        = n_chim_rows,
  n_symbols_pre_hgnc   = n_sym_before,
  n_symbols_post_hgnc  = length(unique(new)),
  biotype_exclude      = params$biotype_exclude,
  n_genes_excluded     = n_excluded,
  excluded_fraction    = list(min = min(excl_frac), median = median(excl_frac),
                              max = max(excl_frac)),
  top_share_retained   = as.list(round(top_share, 5)),
  n_genes_all          = length(symbols_all),
  n_genes_kept         = nrow(tpm),
  n_markers_rescued    = n_rescued,
  min_marker_fpkm      = min_fpkm,
  top_n                = top_n,
  marker_survival      = marker_survival
)
write_json(summary_list, file.path(out_dir, "02_expression_summary.json"),
           auto_unbox = TRUE, pretty = TRUE)

cat("\n== Expression matrix ==\n")
cat(sprintf("samples            : %d\n", ncol(tpm)))
cat(sprintf("transcripts read   : %d\n", n_tx_total))
cat(sprintf("  named directly   : %d (%.1f%%) -- gene_name in the ctab\n",
            n_named, 100*n_named/n_tx_total))
cat(sprintf("  inherited symbol : %d (%.1f%%) -- novel isoform in a clean locus\n",
            n_inherited, 100*n_inherited/n_tx_total))
cat(sprintf("  dropped, no ref  : %d (%.1f%%) -- unnamed, locus has no gene\n",
            n_unnamed, 100*n_unnamed/n_tx_total))
cat(sprintf("  dropped, chimeric: %d (%.1f%%) -- unnamed, locus is fused\n",
            n_chim_rows, 100*n_chim_rows/n_tx_total))
cat(sprintf("symbols            : %d -> %d after HGNC harmonisation\n",
            n_sym_before, length(unique(new))))
if (length(excl_ids)) {
  cat(sprintf("excluded by annotation: %d gene(s) on %s, holding %.1f%%-%.1f%% of FPKM per sample\n",
              n_excluded, paste(params$biotype_exclude, collapse = "/"),
              100 * min(excl_frac), 100 * max(excl_frac)))
} else {
  cat("excluded by annotation: none -- step 1 wrote no exclusion list, so structural\n")
  cat("  and mitochondrial RNA is still in the TPM denominator. Re-run gtf-map with\n")
  cat("  --biotype-tsv if that was not intended.\n")
}
cat(sprintf("genes              : %d all, %d kept by filter (%d markers rescued)\n",
            length(symbols_all), nrow(tpm), n_rescued))
cat(sprintf("  rescue threshold applied on the normalised scale; the same cut on raw\n"))
cat(sprintf("  FPKM would have rescued %d, a difference of %+d marker(s)\n",
            n_rescued_raw, n_rescued - n_rescued_raw))

cat("\n== Largest remaining genes (mean share of retained FPKM) ==\n")
for (i in seq_along(top_share)) {
  cat(sprintf("  %-16s %6.2f%%\n", names(top_share)[i], 100 * top_share[i]))
}
cat("A single gene above ~1% here is worth checking against the annotation: if\n")
cat("its share also varies across samples it distorts the denominator the same way.\n")

cat(sprintf("\n== Marker survival (of top %d per cell type) ==\n", top_n))
print(marker_survival, row.names = FALSE)
cat(sprintf("\nWrote results to %s\n", out_dir))
