#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 4) {
    q()
}

twas.dir <- argv[1]                # twas.dir = './twas/'
depth <- as.integeer(argv[2])      # depth = 7 # tree depth in hierarchical model
prob.cutoff <- as.numeric(argv[3]) # prob.cutoff = 0.05
out.file <- argv[4]                # out.file = 'temp.RData'

library(readr)
library(dplyr)
library(tidyr)
source('Util.R')

.read.tsv <- function(...) suppressMessages(read_tsv(...))

clustering.factors.jaccard <- function(tab, depth = 7, prob.cutoff = .05) {

    ## Construct feature matrix
    feat.mat <- tab %>%
        select(ensg, factor, feature, val) %>%
            spread(key = feature, value = val, fill = 0)

    ## Calculate Jaccard coefficients
    factor.jaccard.mat <- feat.mat %c% (-(1:2)) %>%
        as.matrix() %>%
            (function(x) as(x, 'Matrix::sparseMatrix')) %>%
                text2vec::sim2(method = 'jaccard') %>%
                    as.matrix()

    ## Take cutoff
    jaccard.all.pairs <- which(factor.jaccard.mat > -1, arr.ind = TRUE) %>%
        as.data.frame() %>%
            filter(row > col) %>%
                as.matrix() %>%
                    (function(x) factor.jaccard.mat[x])

    .cutoff <- quantile(jaccard.all.pairs, probs = 1 - prob.cutoff)

    ## Construct network data
    factor.pairs <- which(factor.jaccard.mat > (.cutoff), arr.ind = TRUE) %>%
        as.data.frame() %>%
            filter(row > col)

    edge.weights <- factor.pairs %>%
        as.matrix() %>%
            (function(x) pmin(factor.jaccard.mat[x] / .cutoff, 10)) %>%
                .unlist()

    factor.pairs <- data.frame(factor.pairs, weight = edge.weights)

    V <- nrow(feat.mat)

    hsb.data <- hsblock::pairs.to.sparse.matrix(factor.pairs, vertices = 1:V)
    log.msg('Clustering V = %d factors; E = %.2f edges\n', V, sum(hsb.data$A >= 1)/2)

    ## Run clustering
    hsb <- hsblock::fit.hsblock(hsb.data$A,               # adjacency matrix
                                distrib = 'poisson',      # Poisson model
                                rseed = 20180804,         # some random seed
                                tree.depth = (depth + 1), # need +1 for the root
                                vbiter = 1333,
                                inner.iter = 333,
                                burnin.iter = 10,
                                rate = 0.01,
                                decay = -0.75,
                                delay = 1,
                                record.interval = 3)

    zz <- t(hsb$Z) %>% as.matrix() %>% as.data.frame()
    ret <- data.frame(feat.mat[as.integer(hsb.data$V), 1:2], zz)
    names(ret)[-(1:2)] <- 'K' %&&% 1:ncol(zz)

    list(Z = ret, features = feat.mat, data = hsb.data, cutoff = .cutoff, hsb = hsb)
}

twas.files <- list.files(twas.dir, pattern = 'txt.gz', full.names = TRUE)

twas.tab <- twas.files %>%
    lapply(FUN = .read.tsv) %>%
    bind_rows() %>%
    mutate(trait = gsub(trait, pattern = '-', replacement = '_'))

max.nfactors <- twas.tab %>%
    group_by(trait) %>%
    summarize(nfactors = n()) %>%
    select(nfactors) %>%
    max() %>%
    .unlist()

twas.cutoff <- 0.05 / max.nfactors

selected.factors <- twas.tab %>%
    filter(p.val < twas.cutoff) %>%
    select(ensg, factor) %>%
    unique()

twas.factor.cluster <- twas.tab %>%
    filter(p.val < twas.cutoff) %>%
    mutate(feature = trait, val = 1) %>%
    clustering.factors.jaccard(depth = depth, prob.cutoff = prob.cutoff)

save('twas.factor.cluster', file = out.file)
log.msg('Finished')
