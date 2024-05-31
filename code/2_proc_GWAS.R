# Author: Zheng Li
# Date: 2024-05-30
# Purpose: preprocess GWAS summary statistics from UKBB

library(bigreadr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
gwas_file <- args[1]
chr <- as.numeric(args[2])
setwd("~/VINTAGE_pjt/1_code/VINTAGE-analysis/")
data <- read.table(paste0("data/", gwas_file), header = T, sep = '\t')
snp_info <- do.call(rbind, strsplit(data$variant, split = ":"))
colnames(snp_info) <- c("SNPChr", "SNPPos", "other_allele", "effect_allele")
data <- cbind(data, snp_info)
out <- data %>%
  filter(SNPChr == chr, low_confidence_variant == "false", minor_AF >= 0.01) %>%
  transmute(
    variant = variant,
    chr = SNPChr,
    pos  = SNPPos,
    effect_allele = effect_allele,
    other_allele = other_allele,
    zscore = tstat,
    N = n_complete_samples
  )
fwrite2(out, file = paste0("output/GWAS_chr", chr, "_processed.txt.gz"), 
  col.names = T, row.names = F, sep = " ", compress = "gzip")