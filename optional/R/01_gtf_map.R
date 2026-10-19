#!/usr/bin/env Rscript
##
## 01_gtf_map.R -- locus -> symbol map from the merged StringTie GTF, plus the
## diagnostics deciding whether that assembly can support marker-based
## deconvolution at all.
##
##   Rscript 01_gtf_map.R <params.json>
##
## Exits non-zero only on a hard gate failure (a gated cell type below
## gate_clean_floor), so a SLURM job visibly fails when the answer is
## "re-run StringTie". A soft failure prints a warning and exits 0.
##

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(stringr)
  library(BRETIGEA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript 01_gtf_map.R <params.json>")
params <- fromJSON(args[[1L]])

out_dir    <- params$output_dir
gtf_path   <- params$merged_gtf
top_n      <- params$top_n
gate_min   <- params$gate_clean_min
gate_floor <- params$gate_clean_floor
gate_cts   <- params$gate_cell_types

if (!file.exists(gtf_path)) stop("merged GTF not found: ", gtf_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- transcript lines only -------------------------------------------------
## A human merged GTF is mostly exon lines; filtering in the shell first avoids
## parsing tens of millions of rows we would immediately discard.
reader <- if (grepl("\\.gz$", gtf_path)) "zcat" else "cat"
cmd <- sprintf("%s %s | awk -F'\t' '$0 !~ /^#/ && $3 == \"transcript\"'",
               reader, shQuote(gtf_path))
tx <- fread(cmd = cmd, sep = "\t", header = FALSE, quote = "", showProgress = FALSE)
if (nrow(tx) == 0L) stop("no transcript features found in ", gtf_path)

pull <- function(x, key) str_match(x, paste0(key, ' "([^"]*)"'))[, 2]
txinfo <- data.table(
  locus    = pull(tx$V9, "gene_id"),
  tx_id    = pull(tx$V9, "transcript_id"),
  symbol   = pull(tx$V9, "gene_name"),
  ref_gene = pull(tx$V9, "ref_gene_id")
)
txinfo[symbol %in% c(".", ""), symbol := NA_character_]
if (all(is.na(txinfo$symbol))) {
  stop("no gene_name attributes found in ", gtf_path,
       " -- this looks like a de novo assembly with no reference merged in; ",
       "re-run stringtie -e -B -G <reference.gtf>")
}

## ---- annotation-based exclusion --------------------------------------------
## Genes removed for what they are rather than how they behave, read from the
## list Python derived from --biotype-tsv. This runs BEFORE the chimerism test
## on purpose. A mitochondrial or rRNA gene overlapping a protein-coding one
## makes the merged assembly fuse their loci, and the test below would then
## discard the protein-coding gene -- real, analysable expression lost to a
## collision with something never destined for analysis. Filtering first means
## those loci are never called chimeric.
excl_path <- params$biotype_exclusions
excl_ids  <- character(0)
if (!is.null(excl_path) && nzchar(excl_path) && file.exists(excl_path)) {
  ex <- fread(excl_path, colClasses = "character")
  excl_ids <- unique(c(ex$gene_id, ex$gene_id_base))
  excl_ids <- excl_ids[!is.na(excl_ids) & nzchar(excl_ids)]
}

## The table's identifiers may be ENSG IDs or symbols, versioned or not; all
## three routes are tried rather than requiring one convention of the caller.
txinfo[, ref_base := sub("\\..*$", "", ref_gene)]
is_excluded <- (!is.na(txinfo$ref_gene) & txinfo$ref_gene %in% excl_ids) |
               (!is.na(txinfo$ref_base) & txinfo$ref_base %in% excl_ids) |
               (!is.na(txinfo$symbol)   & txinfo$symbol   %in% excl_ids)

## Measured on the unfiltered assembly, so the report can state what the
## exclusion bought and not only what it cost.
map_pre  <- unique(txinfo[!is.na(symbol), .(locus, symbol)])
chim_pre <- map_pre[, .N, by = locus][N > 1L, locus]
lost_pre <- setdiff(unique(map_pre$symbol),
                    unique(map_pre[!(locus %in% chim_pre), symbol]))

excl_symbols  <- sort(unique(txinfo$symbol[is_excluded & !is.na(txinfo$symbol)]))
n_tx_excluded <- sum(is_excluded)
txinfo <- txinfo[!is_excluded]
if (nrow(txinfo) == 0L) {
  stop("every transcript was excluded by the biotype filter -- check ",
       "--biotype-exclude and the identifier column of --biotype-tsv")
}

## ---- map, and the two failure modes ----------------------------------------
## Chimerism is a property of the LOCUS, not the symbol. Blacklisting every
## symbol that appears at any chimeric locus also discards that gene's clean
## loci elsewhere in the genome, which is real, uncontaminated expression.
## Dropping the loci instead keeps it; only symbols with no clean locus at all
## are genuinely lost.
map        <- unique(txinfo[!is.na(symbol), .(locus, symbol)])
chim_loci  <- map[, .N, by = locus][N > 1L, locus]
clean_map  <- map[!(locus %in% chim_loci)]
resolvable <- unique(clean_map$symbol)
lost       <- sort(setdiff(unique(map$symbol), resolvable))
partial    <- sort(intersect(unique(map[locus %in% chim_loci, symbol]), resolvable))
fragmented <- clean_map[, .N, by = symbol][N > 1L][order(-N)]

## What the exclusion bought: a locus that was fused before and is a clean,
## named locus now, and a symbol that had no clean locus anywhere before and
## does now. Both are protein-coding expression handed back to the analysis.
loci_removed <- setdiff(unique(map_pre$locus), unique(map$locus))
loci_rescued <- setdiff(intersect(chim_pre, unique(map$locus)), chim_loci)
sym_rescued  <- sort(intersect(lost_pre, resolvable))

## ---- marker coverage -------------------------------------------------------
md      <- BRETIGEA::markers_df_human_brain
all_cts <- sort(unique(md$cell))
cts     <- if (is.null(params$cell_types)) all_cts else params$cell_types
unknown <- setdiff(c(cts, gate_cts), all_cts)
if (length(unknown)) {
  stop("cell type(s) not in markers_df_human_brain: ", paste(unknown, collapse = ", "),
       " (available: ", paste(all_cts, collapse = ", "), ")")
}

coverage <- rbindlist(lapply(cts, function(ct) {
  top <- head(md$markers[md$cell == ct], top_n)
  data.table(
    cell_type   = ct,
    n_available = length(top),
    present     = length(intersect(top, unique(map$symbol))),
    clean       = length(intersect(top, resolvable)),
    lost_chim   = length(intersect(top, lost)),
    gated       = ct %in% gate_cts
  )
}))

## The gate keys off `present`, not `clean`: step 2 resolves a transcript's
## symbol from its own gene_name first and only falls back to the locus map for
## novel isoforms, so a marker carrying a gene_name is usable regardless of
## whether its locus is chimeric. `clean` measures something narrower -- which
## markers can additionally inherit novel isoforms -- and is reported, not gated.
gated_cov      <- coverage[gated == TRUE]
gate_pass      <- nrow(gated_cov) > 0L && all(gated_cov$present >= gate_min)
gate_hard_fail <- nrow(gated_cov) > 0L && any(gated_cov$present < gate_floor)

## ---- marker loss by specificity rank ---------------------------------------
## BRETIGEA ranks markers by cell-type specificity, so losing ranks 1-100 costs
## more than losing 400-500. Flat totals cannot tell those cases apart; if
## retention in the top band tracks the bottom band, chimeric loss is rank-blind.
bands <- list(c(1L, 100L), c(101L, 200L), c(201L, 300L), c(301L, 500L))
rank_detail <- rbindlist(lapply(cts, function(ct) {
  ranked <- head(md$markers[md$cell == ct], top_n)
  rbindlist(lapply(bands, function(b) {
    if (b[1] > length(ranked)) return(NULL)
    g <- ranked[seq.int(b[1], min(b[2], length(ranked)))]
    data.table(cell_type = ct,
               rank_band = sprintf("%d-%d", b[1], b[2]),
               n         = length(g),
               clean     = length(intersect(g, resolvable)),
               retention = round(length(intersect(g, resolvable)) / length(g), 3))
  }))
}))

## ---- outputs ---------------------------------------------------------------
## The map written here contains only clean loci, so step 2's join drops
## chimeric loci automatically and needs no separate blacklist.
fwrite(clean_map,  file.path(out_dir, "01_locus_symbol_map.csv"))
fwrite(coverage,   file.path(out_dir, "01_marker_coverage.csv"))
fwrite(rank_detail, file.path(out_dir, "01_marker_rank_bands.csv"))
fwrite(fragmented, file.path(out_dir, "01_fragmented_symbols.csv"))
writeLines(chim_loci, file.path(out_dir, "01_chimeric_loci.txt"))
writeLines(lost,      file.path(out_dir, "01_chimeric_symbols.txt"))
writeLines(partial,   file.path(out_dir, "01_partial_symbols.txt"))
## Step 2 works in symbol space, so the excluded identifiers are translated
## here, where the GTF supplies the ENSG-to-symbol correspondence.
writeLines(excl_symbols, file.path(out_dir, "01_excluded_symbols.txt"))

summary_list <- list(
  merged_gtf           = gtf_path,
  bretigea_version     = as.character(packageVersion("BRETIGEA")),
  r_version            = paste(R.version$major, R.version$minor, sep = "."),
  n_transcripts        = nrow(txinfo),
  n_transcripts_named  = sum(!is.na(txinfo$symbol)),
  biotype_tsv          = params$biotype_tsv,
  biotype_exclude      = params$biotype_exclude,
  n_tx_excluded        = n_tx_excluded,
  n_symbols_excluded   = length(excl_symbols),
  n_chimeric_loci_pre  = length(chim_pre),
  n_loci_rescued       = length(loci_rescued),
  n_loci_removed       = length(loci_removed),
  n_symbols_rescued    = length(sym_rescued),
  n_loci               = uniqueN(txinfo$locus),
  n_loci_named         = uniqueN(map$locus),
  n_loci_clean         = uniqueN(clean_map$locus),
  n_symbols            = uniqueN(map$symbol),
  n_symbols_resolvable = length(resolvable),
  n_chimeric_loci      = length(chim_loci),
  n_chimeric_symbols   = length(lost),
  n_partial_symbols    = length(partial),
  n_fragmented_symbols = nrow(fragmented),
  top_n                = top_n,
  gate_cell_types      = gate_cts,
  gate_clean_min       = gate_min,
  gate_clean_floor     = gate_floor,
  gate_pass            = gate_pass,
  gate_hard_fail       = gate_hard_fail,
  coverage             = coverage,
  rank_detail          = rank_detail
)
write_json(summary_list, file.path(out_dir, "01_gtf_map_summary.json"),
           auto_unbox = TRUE, pretty = TRUE)

## ---- report ----------------------------------------------------------------
cat("\n== GTF map ==\n")
if (length(excl_ids)) {
  cat(sprintf("biotype exclusion  : %d transcript(s), %d symbol(s) removed on %s\n",
              n_tx_excluded, length(excl_symbols),
              paste(params$biotype_exclude, collapse = "/")))
  cat(sprintf("  chimeric loci    : %d before -> %d after, %d locus/loci no longer fused\n",
              length(chim_pre), length(chim_loci), length(loci_rescued)))
  cat(sprintf("  symbols regained : %d had no clean locus before the exclusion and do now\n",
              length(sym_rescued)))
  if (length(sym_rescued)) {
    cat(sprintf("    %s\n", paste(head(sym_rescued, 20L), collapse = ", ")))
    if (length(sym_rescued) > 20L) cat(sprintf("    ... and %d more\n", length(sym_rescued) - 20L))
  }
  cat("Every count below is measured after the exclusion.\n")
} else {
  cat("biotype exclusion  : none (no --biotype-tsv given)\n")
}
cat(sprintf("transcripts        : %d (%d with a gene_name, %.1f%%)\n",
            summary_list$n_transcripts, summary_list$n_transcripts_named,
            100 * summary_list$n_transcripts_named / summary_list$n_transcripts))
cat(sprintf("loci               : %d (%d named, %d clean after dropping chimeric)\n",
            summary_list$n_loci, summary_list$n_loci_named, summary_list$n_loci_clean))
cat(sprintf("symbols            : %d (%d usable from at least one clean locus)\n",
            summary_list$n_symbols, summary_list$n_symbols_resolvable))
cat(sprintf("chimeric loci      : %d -- %d symbol(s) lost entirely, %d kept via a clean locus\n",
            summary_list$n_chimeric_loci, summary_list$n_chimeric_symbols,
            summary_list$n_partial_symbols))
cat(sprintf("fragmented symbols : %d (harmless -- re-joined by symbol in step 2)\n",
            summary_list$n_fragmented_symbols))

cat("\n== Marker coverage (of top", top_n, "per cell type) ==\n")
print(coverage, row.names = FALSE)
cat("`present` is what step 2 can use: a marker with a gene_name resolves\n")
cat("directly, chimeric locus or not. `clean` is the narrower subset that can\n")
cat("also inherit novel isoforms from an unambiguous locus.\n")

cat("\n== Clean markers by specificity rank ==\n")
print(dcast(rank_detail, cell_type ~ rank_band, value.var = "retention"), row.names = FALSE)

cat("\n== Gate ==\n")
if (gate_hard_fail) {
  cat(sprintf("FAIL: a gated cell type has fewer than %d clean markers.\n", gate_floor))
  cat("The merged assembly is too chimeric for marker-based deconvolution.\n")
  cat("Recommended fix: re-run `stringtie -e -B -G <reference.gtf>` on the BAMs.\n")
  cat("Novel transcript discovery is irrelevant here, so nothing is lost.\n")
  quit(status = 1L)
} else if (!gate_pass) {
  cat(sprintf("WARNING: at least one gated cell type is below the %d clean-marker target,\n", gate_min))
  cat(sprintf("but none is below the %d floor. Usable, with the caveat recorded.\n", gate_floor))
} else {
  cat(sprintf("PASS: every gated cell type has at least %d clean markers.\n", gate_min))
}
cat(sprintf("\nWrote results to %s\n", out_dir))
