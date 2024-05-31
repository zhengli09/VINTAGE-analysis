# Author: Zheng Li
# Date: 2024-05-30
# Purpose: preprocess cis-eQTL summary statistics from eQTLGen phase I

library(bigreadr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
eqtl_file <- args[1]
chr <- as.numeric(args[2])
setwd("~/VINTAGE_pjt/1_code/VINTAGE-analysis/")
data <- fread2(paste0("data/", eqtl_file))
out <- data %>% 
  filter(SNPChr == chr) %>%
  transmute(
    gene = Gene,
    variant = paste(SNPChr, SNPPos, OtherAllele, AssessedAllele, sep = ":"),
    chr = SNPChr,
    pos = SNPPos,
    effect_allele = AssessedAllele,
    other_allele = OtherAllele,
    zscore = Zscore,
    N = NrSamples)
fwrite2(out, file = paste0("output/eQTL_chr", chr, "_processed.txt.gz"), 
  col.names = T, row.names = F, sep = " ", compress = "gzip")

