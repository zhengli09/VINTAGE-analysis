# Author: Zheng Li
# Date: 2024-05-30
# Purpose: preprocess gene annotations

library(bigreadr)
library(dplyr)
setwd("~/VINTAGE_pjt/1_code/VINTAGE-analysis/")
args <- commandArgs(trailingOnly = TRUE)
gene_annot <- args[1]
chr <- as.numeric(args[2])

# extract attributes from column 9 of the GENCODE file
get_attrs <- function(x, vars){
  split_data <- lapply(x, function(i){
    sapply(strsplit(i, "\\\";|;")[[1]], function(j){
      elems <- strsplit(j, " \\\"| ")[[1]]
      elems <- elems[elems != ""]
    })
  })
  fields <- lapply(split_data, function(x) x[1, ])
  fields <- unique(unlist(fields))
  
  # initialize an empty data frame with the specified fields
  df <- data.frame(matrix(ncol = length(vars), nrow = length(x)))
  colnames(df) <- vars
  
  # fill the data frame
  for(i in seq_along(split_data)){
    row_data <- split_data[[i]]
    df[i, ] <- row_data[2, match(vars, row_data[1, ])]
  }
  
  return(df)
}

dat <- fread2(paste0("data/", gene_annot), skip = 6)
dat <- subset(dat, V3 == "gene")
gene_attrs <- get_attrs(dat$V9, c("gene_type", "gene_name", "gene_id"))
dat <- cbind(dat, gene_attrs)

out <- dat %>%
  filter(gene_type == "protein_coding") %>%
  transmute(
    gene_id = sapply(strsplit(gene_id, "[.]"), `[`, 1),
    gene_name = gene_name,
    gene_type = gene_type,
    chrom = V1,
    start = V4,
    end = V5,
    strand = V7
  ) %>%
  filter(!chrom %in% c("chrM", "chrX", "chrY")) %>%
  mutate(chrom = substring(chrom, 4)) %>%
  filter(chrom == chr)

fwrite2(out, paste0("output/GENCODE_v40_chr", chr, ".txt"), 
  col.names = T, row.names = F, quote = F, sep = " ")

