#!/usr/bin/env Rscript

source('Util.R')
library(dplyr)
library(readr)

################################################################
## 1. construct null distribution
list.snp.max <- function(...) list.files(..., pattern = '.snp-max.gz', full.names = TRUE) 

read.snp.max.null <- function(dd) {
    .files <- list.snp.max(dd)
    log.msg('# files : %d', length(.files))

    null.fun <- function(x) suppressMessages(read_tsv(x) %>% select(lodds))
    ret <- .files %>%
        lapply(FUN = null.fun) %>%
            .unlist()

    log.msg('Read %s, min %f ~ max %f', dd, min(ret), max(ret))
    return(ret)
}

null.dirs <- 'result/fqtl_null/chr' %&&% 1:22 %&&% '/'

null.stat <- null.dirs %>% lapply(FUN = read.snp.max.null) %>%
    .unlist()


read.snp.max <- function(dd) {
    .files <- list.snp.max(dd)
    log.msg('# files : %d', length(.files))

    .fun <- function(x) suppressMessages(read_tsv(x) %>% select(gene, factor, lodds))
    ret <- .files %>%
        lapply(FUN = .fun) %>%
            bind_rows()

    log.msg('Read %s, min %f ~ max %f', dd, min(ret$lodds), max(ret$lodds))
    return(ret)
}

obs.dirs <- 'result/fqtl_obs/chr' %&&% 1:22 %&&% '/'

obs.stat <- obs.dirs %>% lapply(FUN = read.snp.max) %>%
    bind_rows()

################################################################
## 2. show empirical p-value
gene.cols <- c('chr', 'lb', 'ub', 'strand', 'ensg', 'hgnc')
gene.cols.types <- 'ciiccc_'
gene.file <- 'data/coding.genes.txt.gz'

gene.info <- read_tsv(gene.file, col_names = gene.cols, col_types = gene.cols.types) %>%
    mutate(tss = if_else(strand == '+', lb, ub)) %>%
        mutate(gene = 1:n())

library(qvalue)

null.min <- quantile(obs.stat$lodds, 0.05)
null.stat.boot <- sample(pmax(null.stat, null.min), size = nrow(obs.stat) * 5, replace = TRUE)
pval <- qvalue::empPvals(obs.stat$lodds, null.stat.boot)
qval <- qvalue(pval)
qval.bh <- p.adjust(pval, method = 'fdr')

obs.stat.pq <- obs.stat %>% mutate(pval, qval.emp = qval$qvalues, qval.bh)

take.cutoff <- function(tab) {
    cutoff <- tab$lodds
    obs.stat.pq %>% filter(lodds >= cutoff) %>%
        arrange(desc(pval)) %>% head(1) %>%
            mutate(lodds) %>%
                select(lodds, pval, qval.emp, qval.bh)
}

fdr.tab <- data.frame(pip = c(1e-4, 1e-3, 1e-2, seq(.1, .9, .1))) %>%
    mutate(lodds = log(pip) - log(1 - pip)) %>%
        group_by(pip) %>%
            do(take.cutoff(.))

dir.create('stat')

write_tsv(fdr.tab, path = 'stat/fdr.txt.gz')

################################################################
## 3. combine all the results chromosome by chromosome
## with logit pip cutoff > -2.18

list.combined <- function(...) list.files(..., pattern = '.combined.gz', full.names = TRUE) 

valid.genes <- obs.stat %>% filter(lodds > -2.18) %>%
    select(gene, factor) %>% unique()

read.dd <- function(dd) {
    ret <- list.combined(dd) %>%
        lapply(FUN = function(x) suppressMessages(read_tsv(x, col_types = 'iicccccccc'))) %>%
            bind_rows()
    log.msg('Read %d rows : %s', nrow(ret), dd)
    return(ret)
}

dir.create('stat')

for(chr in 1:22) {
    out.tab <- read.dd(obs.dirs[chr])
    ret <- gene.info %>%
        left_join(valid.genes) %>%
            left_join(out.tab) %>%
                na.omit() %>%
                    select(-gene)                

    write_tsv(ret, path = 'stat/fqtl_' %&&% chr %&&% '.txt.gz')
}

