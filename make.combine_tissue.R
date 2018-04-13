#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
source('Util.FQTL.R')
library(readr)
library(dplyr)
library(tidyr)

if(length(argv) != 2) q()

job.id <- as.integer(argv[1])
out.dir <- argv[2]
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

################################################################
## organize tissue samples

tis.sample.info <- read.tis.sample.info()

tis.tab <- tis.sample.info %>%
    select(tis.idx, tis.dir, tis.name) %>%
        unique() %>%
            arrange(tis.idx)


################################################################
## write out gene-by-gene
data.dir <- 'processed/expression/' %&&% tis.tab$tis.dir
data.files <- data.dir %&&% ('/nb_' %&&% job.id %&&% '.txt.gz')


.read.tis <- function(j) {
    ret <- suppressMessages(read_tsv(data.files[j]))
    if(nrow(ret) == 0) return(NULL)
    ret <- ret %>%
        mutate(expr.pos = 1:n()) %>%
            gather(key = gene.idx, value = expr, -expr.pos) %>%
                mutate(tis.idx = j, gene.idx = as.integer(gene.idx))
    return(ret)
}

expr.data <- lapply(1:length(data.files), FUN = .read.tis) %>%
    bind_rows() %>%
        left_join(tis.sample.info %>% select(tis.idx, expr.pos, geno.pos)) %>%
            na.omit() %>%
                select(-expr.pos) %>%
                    as.data.frame()

genes <- expr.data$gene.idx %>% unique()

## construct individual x tissue matrix and store them
for(gg in genes) {
    temp <- expr.data %>% mutate(expr = round(expr, 3)) %>%
        filter(gene.idx == gg) %>%
            spread(key = tis.idx, value = expr)

    out.file <- out.dir %&&% '/' %&&% gg %&&% '.txt.gz'
    write_tsv(temp, path = out.file, col_names = TRUE)
    log.msg('Wrote %s\n', out.file)
    rm(temp)
}

log.msg('Successfully completed\n')
