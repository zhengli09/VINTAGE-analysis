# Author: Zheng Li
# Date: 2024-05-30
# Purpose: preprocess LD reference data from 1000G

setwd("~/VINTAGE_pjt/1_code/VINTAGE-analysis/")
ld_ref_file <- commandArgs(trailingOnly = TRUE)
ld_ref_file <- paste0("data/", ld_ref_file)

# find European sample IDs
meta_file <- "data/integrated_call_samples_v3.20130502.ALL.panel"
meta <- read.table(meta_file, header = T)
eur_smps <- meta$sample[meta$super_pop == "EUR"]
smps_file <- "output/1kg_eur_smps.txt"
write.table(eur_smps, file = smps_file, col.names = F, row.names = F, quote = F)

# extract genotype data of the European samples
ld_out_file <- "output/LD_1kg_chr22_processed"
cmd <- paste("/usr/cluster/bin/plink2 --vcf", ld_ref_file, "--keep", smps_file, 
  "--snps-only --maf 0.05 --max-alleles 2 --make-bed --out", ld_out_file)
system(cmd)
