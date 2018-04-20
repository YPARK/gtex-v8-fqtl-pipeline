#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
library(readr)
library(dplyr)
library(tidyr)
library(fqtl)

if(length(argv) != 4) { q() }

expr.file <- argv[1]                  # e.g., expr.file = 'processed/combined/chr1/4.txt.gz'
geno.hdr <- argv[2]                   # e.g., geno.hdr = 'geno/chr1'
out.hdr <- argv[3]                    # e.g., out.hdr = 'temp'
do.permutation <- as.logical(argv[4]) # e.g., TRUE

################################################################
temp.dir0 <- '/broad/hptmp/ypp/gtex/v8/' %&&% out.hdr

cmd <- 'mkdir -p ' %&&% temp.dir0 %&&% '; mktemp -d ' %&&% temp.dir0 %&&% '_temp.XXXXXX'
temp.dir <- system(cmd, intern = TRUE)

cis.dist <- 1e6
rseed <- 1667

################################################################

y.tab <- read_tsv(expr.file)
gene.idx <- unique(y.tab$gene.idx)
stopifnot(length(gene.idx) == 1)

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

joint.out.file <- out.hdr %&&% '.combined.gz'
snp.out.file <- out.hdr %&&% '.snp-factor.gz'
max.snp.out.file <- out.hdr %&&% '.snp-max.gz'
tis.out.file <- out.hdr %&&% '.tis-factor.gz'

if(file.exists(joint.out.file)) {
    log.msg('File already exists: %s', joint.out.file)
    q()
}

################################################################
## 1. Find cis genotype matrix
gene.cols <- c('chr', 'lb', 'ub', 'strand', 'ensg', 'hgnc', 'remove')
gene.file <- 'data/coding.genes.txt.gz'

gene.info <- read_tsv(gene.file, col_names = gene.cols) %>%
    dplyr::select(-remove) %r% gene.idx %>%
        mutate(tss = if_else(strand == '+', lb, ub))

chr <- gene.info$chr %>% .unlist()
lb <- max(gene.info$tss - cis.dist, 0)
ub <- max(gene.info$tss + cis.dist, 0)

dir.create(temp.dir, recursive = TRUE)
plink <- subset.plink(geno.hdr, chr, lb, ub, temp.dir)
system('rm -r ' %&&% temp.dir)

################################################################
## 2. fit the FQTL model

geno.pos <- y.tab$geno.pos
X <- plink$BED %r% geno.pos %>% scale() %>% rm.na.zero() %>% as.matrix()
Y <- y.tab %>% select(- gene.idx, -geno.pos) %>% scale() %>% as.matrix()
tis.idx <- colnames(Y)

if(do.permutation) {
    set.seed(rseed)
    n <- dim(Y)[1]

    Y <- apply(Y, 2,
               function(y) {
                   ret <- matrix(NA, n, 1)
                   ret.pos <- is.finite(y)
                   y.shuf <- y[ret.pos]
                   y.shuf <- sample(y.shuf)
                   ret[ret.pos, 1] <- y.shuf
                   return(ret)
               })

    log.msg('Permuted breaking tissue-tissue correlation\n')
}

K <- min(ncol(Y), ncol(X))

opt.reg <- list(vbiter = 5000, gammax = 1e4, tol = 1e-8,
                rate = 1e-2, decay = -1e-2,
                pi.ub = -1/2, pi.lb = -2, tau = -4, do.hyper = TRUE,
                jitter = 0.1, svd.init = TRUE, out.residual = FALSE,
                print.interv = 10, k = K)

fqtl.out <- fqtl.regress(y = Y, x.mean = X, factored = TRUE, options = opt.reg)

log.msg('Successfully finished model estimation')

snp.lodds.cutoff <- log(0.1) - log(0.9)

tis.effect <- fqtl.out$mean.right %>% melt.spike.slab() %>%
    rename(tis = row, factor = col) %>%
        mutate(gene = gene.idx) %>%
            select(gene, tis, factor, lodds, theta, theta.sd)

bim.tab <- plink$BIM %>% select(rs) %>% mutate(snp = 1:n())

snp.effect <- fqtl.out$mean.left %>% melt.spike.slab() %>%
    rename(snp = row, factor = col) %>%
        mutate(gene = gene.idx) %>%
            left_join(bim.tab) %>%
                select(gene, rs, factor, lodds, theta, theta.sd)

snp.max.effect <- snp.effect %>%
    group_by(gene, factor) %>%
        slice(which.max(lodds))

snp.effect <- snp.effect %>%
        filter(lodds > snp.lodds.cutoff)

out <- tibble()

.collapse <- function(...) paste(..., collapse = '|')

if(nrow(snp.effect) > 0) {

    out <- snp.effect %>% group_by(gene, factor) %>%
        summarize(snp = .collapse(rs),
                  snp.theta = .collapse(theta),
                  snp.sd = .collapse(theta.sd),
                  snp.lodds = .collapse(lodds))

    temp <- tis.effect %>% group_by(gene, factor) %>%
        summarize(tis = .collapse(tis),
                  tis.theta = .collapse(theta),
                  tis.sd = .collapse(theta.sd),
                  tis.lodds = .collapse(lodds))

    out <- out %>% left_join(temp)
}

write_tsv(tis.effect, path = tis.out.file)
write_tsv(snp.effect, path = snp.out.file)
write_tsv(out, path = joint.out.file)
write_tsv(snp.max.effect, path = max.snp.out.file)

log.msg('Successfully finished everything!')
