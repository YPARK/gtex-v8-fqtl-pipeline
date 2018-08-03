#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

twas.dir <- argv[1] # 'twas/imputed_ADIPOGen_Adiponectin'
out.file <- argv[2]

options(stringsAsFactors = FALSE)
source('Util.R')
library(dplyr)
library(tidyr)
library(readr)

trait.name <- basename(twas.dir) %>%
    gsub(pattern = 'imputed_', replacement = '')

list.twas.files <- function(x) {
    list.files(path = x, pattern = '^chr.*\\.txt\\.gz$', full.names = TRUE)
}
.read.tsv <- function(x) {
    ret <- suppressMessages(read_tsv(x))
    chr <- basename(x) %>% strsplit(split = '_') %>%
        (function(s) s[[1]][1])
    if(nrow(ret) < 1) return(NULL)
    return(ret %>% mutate(chr))
}

coding.cols <- c('chr', 'lb', 'ub', 'strand', 'ensg', 'hgnc')
coding.genes <- read_tsv('data/coding.genes.txt.gz',
                         col_names = coding.cols,
                         col_types = 'ciiccc_')

twas.tab <- twas.dir %>% list.twas.files() %>%
    lapply(FUN=.read.tsv) %>% bind_rows() %>%
        arrange(p.val, desc(abs(z))) %>%
            left_join(coding.genes) %>%
                mutate(trait = trait.name)

write_tsv(twas.tab, out.file)
