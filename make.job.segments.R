#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) {
    q()
}


gene.file <- argv[1] # e.g., gene.file = 'data/coding.genes.txt.gz'
out.file <- argv[2]  # e.g., 'jobs.txt'

options(stringsAsFactors = FALSE)
source('Util.R')
library(dplyr)
library(readr)

gene.cols <- c('chr', 'lb', 'ub', 'strand', 'ensg', 'hgnc')
gene.types <- 'ciiccc_'

## 1. assign genes into blocks -- never change the order
genes <-
    read_tsv(gene.file, col_names = gene.cols, col_types = gene.types) %>%
        mutate(loc = if_else(strand == '+', lb, ub)) %>%
            mutate(data.pos = 1:n())

block.size0 <- 1e8
max.size <- 50

genes.blk <-
    genes %>%
        mutate(block = floor(loc/block.size0), depth = 0) %>%
            select(chr, loc, block, data.pos, depth)

big.blocks <- genes.blk %>%
    group_by(chr, block) %>%
        summarize(n.genes = n()) %>%
            filter(n.genes > max.size) %>%
                select(chr, block, n.genes)

## 2. break large blocks into smaller chunks
while(TRUE) {

    big.blocks <- genes.blk %>%
        group_by(chr, block) %>%
            summarize(n.genes = n()) %>%
                filter(n.genes > max.size) %>%
                    select(chr, block)

    if(nrow(big.blocks) == 0) break

    ret1 <- genes.blk %>%
        anti_join(big.blocks)

    ret2 <-
        genes.blk %>% right_join(big.blocks) %>%
            mutate(depth = depth + 1) %>%
                mutate(block = floor(loc / block.size0 * 2^depth)/2^depth)

    genes.blk <- bind_rows(ret1, ret2) %>%
        arrange(data.pos)
}

ret <- genes %>%
    left_join(genes.blk) %>%
        group_by(chr, block) %>%
            summarize(lb = min(data.pos), ub = max(data.pos)) %>%
                arrange(lb) %>%
                        select(chr, lb, ub)

write.table(ret, file = out.file, quote = FALSE, sep = '\t',
            col.names = FALSE, row.names = FALSE)
