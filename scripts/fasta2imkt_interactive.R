#!/usr/bin/env Rscript
# ------------------------------------------------------------
# FASTA → iMKT .daf and .div generator (0.025 bins, reference+ingroup+outgroup)
# USAGE: Rscript fasta2imkt.R input.fasta output_prefix
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Biostrings)
})

if (interactive()) {
  # Please edit this line
  args <- c("./data/jingwei_dteissieri_dyakuba/jingwei_dteissieri_dyakuba_aligned.fasta", "jingwei_dteissieri_dyakuba_aligned")
} else {
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args) < 2) {
  cat("Usage: Rscript fasta2imkt.R input.fasta output_prefix\n")
  quit(status = 1)
}

fasta_file <- args[1]
output_prefix <- args[2]

# ===== Load sequences =====
seqs <- readDNAStringSet(fasta_file)
if (length(seqs) < 3) stop("Need ≥3 sequences (reference, ≥1 ingroup, outgroup).")

ref <- toupper(as.character(seqs[[1]]))
out <- toupper(as.character(seqs[[length(seqs)]]))
ingroup <- toupper(as.character(seqs[2:(length(seqs)-1)]))

# ===== Settings =====
ref <- substr(ref, 1, nchar(ref) - nchar(ref) %% 3)
out <- substr(out, 1, nchar(out) - nchar(out) %% 3)
ingroup <- lapply(ingroup, function(s) substr(s, 1, nchar(ref)))

ncodons <- nchar(ref) / 3
bins <- seq(0.025, 0.975, by = 0.05)

# ===== Codon translation setup =====
std_table <- GENETIC_CODE
stop_codons <- names(std_table)[std_table == "*"]
aa_table <- std_table[std_table != "*"]

translate_codon <- function(codon) {
  codon <- toupper(codon)
  if (grepl("N", codon) || codon %in% stop_codons) return(NA)
  aa <- aa_table[[codon]]
  if (is.null(aa)) return(NA)
  return(aa)
}

compare_codons <- function(a, b) {
  if (a == b || grepl("N", paste0(a, b))) return(c(0, 0))
  aa_a <- translate_codon(a)
  aa_b <- translate_codon(b)
  if (is.na(aa_a) || is.na(aa_b)) return(c(0, 0))
  if (aa_a == aa_b) return(c(1, 0)) else return(c(0, 1))
}

# ===== .div calculation =====
syn_div <- 0
nonsyn_div <- 0
for (i in seq_len(ncodons)) {
  cod_ref <- substr(ref, (i - 1) * 3 + 1, i * 3)
  cod_out <- substr(out, (i - 1) * 3 + 1, i * 3)
  res <- compare_codons(cod_ref, cod_out)
  syn_div <- syn_div + res[1]
  nonsyn_div <- nonsyn_div + res[2]
}

# ===== .daf binning =====
Pi <- numeric(length(bins))
P0 <- numeric(length(bins))

for (i in seq_len(ncodons)) {
  cod_ref <- substr(ref, (i - 1) * 3 + 1, i * 3)
  codon_all <- c(cod_ref, sapply(ingroup, function(seq) substr(seq, (i - 1) * 3 + 1, i * 3)))
  unique_codons <- unique(codon_all)
  if (length(unique_codons) == 1) next
  aas <- unique(na.omit(sapply(unique_codons, translate_codon)))
  freq <- sum(codon_all != cod_ref) / length(codon_all)
  idx <- which.min(abs(bins - freq))
  if (length(aas) > 1) Pi[idx] <- Pi[idx] + 1 else P0[idx] <- P0[idx] + 1
}

# ===== Write DAF =====
daf_path <- paste0(output_prefix, ".daf")
daf_df <- data.frame(daf = sprintf("%.3f", bins), Pi = Pi, P0 = P0)
write.table(daf_df, daf_path, quote = FALSE, sep = "\t", row.names = FALSE)

# ===== Write DIV =====
div_path <- paste0(output_prefix, ".div")
div_df <- data.frame(mi = ncodons, Di = nonsyn_div, m0 = ncodons, D0 = syn_div)
write.table(div_df, div_path, quote = FALSE, sep = "\t", row.names = FALSE)

cat(paste0("✔ Created ", daf_path, "\n"))
cat(paste0("✔ Created ", div_path, "\n"))

