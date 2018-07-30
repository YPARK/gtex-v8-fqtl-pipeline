#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 5) {
    q()
}

fqtl.stat.file <- argv[1]             # e.g., fqtl.stat.file = 'stat/fqtl_22.txt.gz'
gwas.file <- argv[2]                  # e.g., gwas.file = 'gwas_data/ADIPOGen_Adiponectin/22.txt.gz'
geno.hdr <- argv[3]                   # e.g., geno.hdr = 'geno/chr22'
gg.vec <- eval(parse(text = argv[4])) # e.g., gg.vec = 1:5
out.file <- argv[5]                   # e.g., out.file = 'temp.txt.gz'

################################################################
source('Util.R')
source('Util.TWAS.R')
library(readr)
library(dplyr)
library(tidyr)
library(zqtl)

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

if(file.exists(out.file)) {
    log.msg('File exists: %s', out.file)
    q()
}

if(!file.exists(gwas.file)) {
    log.msg('GWAS file does not exist: %s', gwas.file)
    q()
}

if(!file.exists(fqtl.stat.file)) {
    log.msg('stat file does not exist: %s', fqtl.stat.file)
    q()
}

################################################################

cis.dist <- 1e6
n.cutoff <- 10
n.perm <- 5e7
n.blk <- 2048
n.round <- ceiling(n.perm/n.blk)

temp.dir <- system('mkdir -p /broad/hptmp/ypp/gtex-v8-twas/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/gtex-v8-twas/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

if(!dir.exists(temp.dir)) {
    log.msg('Failed to create a temporary directory: %s', temp.dir)
    q()
}

fqtl.stat <- read_tsv(fqtl.stat.file)

snp.stat <- fqtl.stat %>%
    select(ensg, factor, snp, snp.theta, snp.sd, snp.lodds)

gene.loc <- fqtl.stat %r% gg.vec %>% select(chr, tss) %>%
    rename(chromosome = chr)

gwas.tab <- read_tsv(gwas.file) %>%    
    filter(chromosome %in% unique(gene.loc$chromosome),
           position >= (min(gene.loc$tss) - cis.dist),
           position <= (max(gene.loc$tss) + cis.dist))
gc()

if(nrow(gwas.tab) < 1) {
    write_tsv(data.frame(), path = out.file)
    log.msg('Finished: empty GWAS')
    q()
}

twas.factor <- function(gg, lodds.cutoff = log(0.9) - log(0.1)) {

    tss <- fqtl.stat[gg, 'tss'] %>% .unlist() %>% as.integer()
    info <- fqtl.stat %r% gg %>% select(ensg, factor)

    plink <- subset.plink(geno.hdr, chr = basename(geno.hdr),
                          max(tss - cis.dist, 0),
                          tss + cis.dist, temp.dir)

    if(is.null(plink)) return(NULL)

    .split <- function(...) .unlist(...) %>% strsplit(split = '[|]') %>% (function(x) x[[1]])

    eqtl.tab <-
        data.frame(snp = snp.stat[gg, 'snp'] %>% .split(),
                   eqtl.beta = snp.stat[gg, 'snp.theta'] %>% .split() %>% as.numeric(),
                   eqtl.se = snp.stat[gg, 'snp.sd'] %>% .split() %>% as.numeric(),
                   eqtl.lodds = snp.stat[gg, 'snp.lodds'] %>% .split() %>% as.numeric())

    eqtl.tab <- eqtl.tab %>% filter(eqtl.lodds >= lodds.cutoff) %>%
        separate(snp, c('chromosome', 'snp.loc', 'eqtl.a1', 'eqtl.a2'), sep = ':') %>%
            mutate(snp.loc = as.integer(snp.loc))
    
    matched <- match.eqtl.gwas(plink, gwas.tab, eqtl.tab)

    if(nrow(matched) < 1) return(NULL)

    eqtl.z.tab <- matched %>% filter(!is.na(eqtl.z))

    if(nrow(eqtl.z.tab) < 1) return(NULL)

    svd.out <- zqtl::take.ld.svd(plink$BED, options = list(eig.tol = 1e-2, do.stdize = TRUE))

    V.t <- svd.out$V.t %c% matched$x.pos
    D <- svd.out$D

    V.t.obs <- svd.out$V.t %c% eqtl.z.tab$x.pos
    gwas.z <- eqtl.z.tab$gwas.z
    eqtl.z <- eqtl.z.tab$eqtl.z
    gwas.z.tot <- matched$gwas.z

    obs.stat <- func.NWAS(eqtl.z, gwas.z, V.t.obs, D)

    blk.mat <- func.blk.ind(nrow(eqtl.z.tab), n.blk)

    ## adaptive permutation
    c.tot <- 0
    p.tot <- 0
    z.obs.abs <- abs(obs.stat[, 'z'])
    for(b in seq(1, n.round)){
        set.seed(b)

        stat.z <- func.NWAS.eqtl.perm(eqtl.z, gwas.z.tot, V.t, D, blk.mat) %>%
            select(z) %>% .unlist()

        c.tot <- c.tot + sum((abs(stat.z) + 1e-8) >= z.obs.abs)
        p.tot <- p.tot + length(stat.z)

        if(c.tot > n.cutoff){
            break;
        }
        cat('\n', c.tot, '/', p.tot, '... ')
    }

    log.msg('Finished: %d, %s', gg, info[1])

    ret <- bind_cols(info, obs.stat) %>% mutate(p.val = (1 + c.tot)/(1 + p.tot))
    return(ret)
}

out.tab <- lapply(gg.vec, twas.factor) %>% bind_rows()

if(nrow(out.tab) < 1) {
    if(dir.exists(temp.dir)) { system('rm -r ' %&&% temp.dir) }
    write_tsv(data.frame(), path = out.file)
    log.msg('Finished: empty results')
    q()
}

out.tab <- out.tab %>%
    mutate(z = signif(z, 4),
           theta = signif(theta, 4),
           theta.se = signif(theta.se, 4),
           p.val = signif(p.val, 4))

if(dir.exists(temp.dir)) { system('rm -r ' %&&% temp.dir) }
write_tsv(out.tab, path = out.file)
log.msg('Finished')
