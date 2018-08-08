#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 5) {
    q()
}

fqtl.stat.file <- argv[1]             # e.g., fqtl.stat.file = 'stat/fqtl_22.txt.gz'
gwas.file <- argv[2]                  # e.g., gwas.file = 'gwas_data/imputed_ADIPOGen_Adiponectin/22.txt.gz'
geno.hdr <- argv[3]                   # e.g., geno.hdr = 'geno/chr22'
gg.vec <- eval(parse(text = argv[4])) # e.g., gg.vec = 1:5
out.file <- argv[5]                   # e.g., out.file = 'temp.txt.gz'

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

cis.dist <- 1e6

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

twas.data.factor <- function(gg, lodds.cutoff) {

    tss <- fqtl.stat[gg, 'tss'] %>% .unlist() %>% as.integer()
    info <- fqtl.stat %r% gg %>% select(hgnc, ensg, factor)
    
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

    gwas.info <- matched %>%
        slice(which.min(gwas.p)) %>%
            select(snp.loc, gwas.p)

    svd.out <- zqtl::take.ld.svd(plink$BED, options = list(eig.tol = 1e-2, do.stdize = TRUE))
    V.t <- svd.out$V.t %c% eqtl.z.tab$x.pos
    D <- svd.out$D
    gwas.z <- eqtl.z.tab$gwas.z
    eqtl.z <- eqtl.z.tab$eqtl.z
    obs.stat <- func.NWAS(eqtl.z, gwas.z, V.t, D)

    out <- eqtl.z.tab %>%
        select(-x.pos) %>%
            mutate(hgnc = .unlist(info$hgnc),
                   ensg = .unlist(info$ensg),
                   factor = as.integer(info$factor),
                   twas.z = signif(as.numeric(obs.stat$z), 4),
                   best.gwas.loc = as.integer(gwas.info$snp.loc),
                   best.gwas.p = signif(as.numeric(gwas.info$gwas.p), 4))

}

out.tab <- lapply(gg.vec, twas.data.factor, lodds.cutoff = 0) %>%
    bind_rows()

if(nrow(out.tab) < 1) {
    if(dir.exists(temp.dir)) { system('rm -r ' %&&% temp.dir) }
    write_tsv(data.frame(), path = out.file)
    log.msg('Finished: empty results')
    q()
}

if(dir.exists(temp.dir)) { system('rm -r ' %&&% temp.dir) }
write_tsv(out.tab, path = out.file)
log.msg('Finished')
