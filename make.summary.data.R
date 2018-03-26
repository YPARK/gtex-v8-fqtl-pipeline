#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) {
    q()
}

options(stringsAsFactors = FALSE)
source('Util.R')
library(fqtl)
library(dplyr)
library(readr)
library(methods)

chr <- as.integer(argv[1])             # e.g., chr = 19
gene.idx.str <- argv[2]                # e.g., gene.idx.str = '1:33'
tis.name <- argv[3]                    # e.g., tis.name = '10_Brain_Cerebellar_Hemisphere'
n.ctrl.per.gene <- as.integer(argv[4]) # e.g., n.ctrl.per.gene = 3
data.type <- argv[5]                   # e.g., data.type = 'rsem-count'
out.hdr <- argv[6]                     # e.g., out.hdr = 'temp'

out.resid.file <- glue(out.hdr, '.resid.gz')
out.qtl.file <- glue(out.hdr, '.qtl.gz')
out.sample.file <- glue(out.hdr, '.samples.gz')
out.gene.file <- glue(out.hdr, '.genes.gz')

gene.idx <- eval(parse(text = gene.idx.str))

################################################################
## scratch pad
temp.dir <- system(glue('mkdir -p ', out.hdr, '; mktemp -d ', out.hdr, '/temp.XXXXXXXX'),
                   intern = TRUE, ignore.stderr = TRUE)


################################################################
## take expression matrices
data.file <- glue('data/', tis.name, '/chr', chr, '/', data.type, '.txt.gz')
gene.file <- glue('data/', tis.name, '/chr', chr, '/genes.txt.gz')
sample.file <- glue('data/', tis.name, '/samples.txt.gz')

other.data.files <- glue('data/', tis.name, '/chr', setdiff(1:22, chr), '/rsem-count.txt.gz')
size.factor.file <- glue('data/', tis.name, '/rsem-size-factor.txt.gz')

.read.tsv <- function(...) { readr::read_tsv(..., col_names = FALSE) %>% as.matrix() }

################################################################
## transformation of raw data into 
size.factor <- .read.tsv(size.factor.file)
adjust.size <- function(mat, sf = size.factor[, 1]) { sweep(mat, 1, sf, `/`) }

rm.zero <- function(mat) { mat[mat == 0] <- NA; return(mat) }

stdize.count <- function(xx) {
    xx[xx < .5] <- NA
    xx.med <- apply(xx, 2, median, na.rm = TRUE)
    xx.scaled <- sweep(xx, 2, xx.med, `/`)
    ret <- xx.scaled * 10
}

.trans.normal <- function(..., .data.type = c('rsem-count', 'vsn')) {
    .data.type <- match.arg(.data.type)
    if(.data.type == 'rsem-count') {
        ret <- log2(1/2 + rm.zero(...) %>% adjust.size()) %>% scale()        
    } else if(.data.type == 'vsn') {
        ret <- scale(...)
    }
    return(ret)
}

.pre <- function(..., .data.type = c('rsem-count', 'vsn')) {
    .data.type <- match.arg(.data.type)
    if(.data.type == 'rsem-count') {
        ret <- rm.zero(...) %>% adjust.size() %>% stdize.count()
    } else if(.data.type == 'vsn') {
        ret <- scale(...)
    }
    return(ret)    
}

.model.name <- function(.data.type = c('rsem-count', 'vsn')) {
    .data.type <- match.arg(.data.type)
    if(.data.type == 'rsem-count') {
        ret <- 'nb'
    } else if(.data.type == 'vsn') {
        ret <- 'gaussian'
    }
    return(ret)        
}


trans.normal <- function(...) .trans.normal(..., .data.type = data.type)
model.name <- .model.name(data.type)
preprocess <- function(...) .pre(..., .data.type = data.type)

find.most.cor.genes <- function(.others, Y, n.ctrl = n.ctrl.per.gene, na.permit = 0) {
    valid.others <- which(apply(is.na(.others), 2, sum) <= na.permit)
    .others <- .others %c% valid.others
    .cor <- abs(fast.cor(.others, Y))
    .ret <- lapply(1:ncol(.cor), function(j) order(.cor[, j], decreasing = TRUE)[1:n.ctrl])
    .ret <- valid.others[unique(do.call(c, .ret))]
    return(.ret)
}

take.cor.genes.chr <- function(other.chr, Y, n.ctrl = n.ctrl.per.gene) {
    .others <- .read.tsv(other.data.files[other.chr])
    .ret <- find.most.cor.genes(.others %>% trans.normal(), Y %>% trans.normal(), n.ctrl)
    return(.others %c% .ret)
}

subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {
    require(fqtl)
    require(dplyr)

    plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, chr, plink.lb, plink.ub, glue(temp.dir, '/plink'))
    system(plink.cmd)

    plink <- read.plink(glue(temp.dir, '/plink'))
    colnames(plink$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
    colnames(plink$FAM) <- c('fam', 'iid', 'father', 'mother', 'sex.code', '.pheno')
    plink$FAM <- plink$FAM %>% mutate(iid = sapply(iid, gsub, pattern = 'GTEX-', replacement = ''))
    return(plink)
}

################################################################
gene.cols <- c('chr', 'tss', 'tes', 'strand', 'ensg', 'hgnc', 'remove')
genes <- read.table(gene.file, col.names = gene.cols) %>%
    dplyr::select(-remove) %r% gene.idx

sample.cols<-c('sampid','iid','data.pos','RIN','TSISCH','tis.name',
               'SEX','AGE','HGHT','WGHT','BMI','SEX.std','AGE.std','HGHT.std',
               'WGHT.std','BMI.std','ntis','tis.idx','tis.dir')

samples <- read.table(sample.file, col.names = sample.cols) %>%
    select(-data.pos, -tis.name, -tis.idx, -tis.dir,
           -SEX, -AGE, -HGHT, -WGHT, -BMI, -ntis, -RIN, -TSISCH)
           
cis.dist <- 1e6
lb <- pmax(min(genes$tss) - cis.dist, 0)
ub <- pmax(max(genes$tes) + cis.dist, 0)
plink <- subset.plink(glue('geno/chr', chr), chr, lb, ub, temp.dir) 

x.fam <- plink$FAM %>% select(iid) %>% mutate(x.pos = 1:n())
x.bim <- plink$BIM

samples.matched <- samples %>% mutate(y.pos = 1:n()) %>%
    left_join(x.fam, by = 'iid') %>% na.omit()

x.pos <- samples.matched$x.pos
y.pos <- samples.matched$y.pos

################################################################
Y <- .read.tsv(data.file) %c% gene.idx
Y.conf <- do.call(cbind, lapply(1:length(other.data.files), take.cor.genes.chr, Y = Y))
gc()

.conf <- find.most.cor.genes(Y.conf %>% trans.normal(), Y %>% trans.normal(), n.ctrl.per.gene)
Y.conf <- Y.conf %c% .conf

################################################################

Y.std <- Y %>% preprocess() %r% y.pos
Y.conf.std <- Y.conf %>% trans.normal() %r% y.pos
X.std <- plink$BED %r% x.pos %>% scale()

C <- samples.matched %>%
    dplyr::select(contains('std')) %>%
        mutate(intercept = 1) %>%
        as.matrix()

vb.opt <- list(vbiter = 5000,
               gammax = 1e4,
               out.residual = TRUE,
               tol = 1e-8,
               rate = 1e-2,
               decay = -1e-2,
               pi = -2,
               tau = -4,
               do.hyper = FALSE,
               model = model.name,
               svd.init = TRUE,
               jitter = 0.01,
               k = 5)

fqtl.out <- fqtl.regress(y = Y.std, x.mean = Y.conf.std, c.mean = C,
                         factored = TRUE, options = vb.opt)
                         
## generate eQTL summary statistics
R <- fqtl.out$resid$theta
R[is.na(Y.std)] <- NA
R <- R %>% scale()
colnames(X.std) <- x.bim$rs
colnames(R) <- genes$ensg
rownames(R) <- samples.matched$iid

rm.ensg <- genes$ensg[apply(is.na(R), 2, mean) > .5]

marginal.qtl <- get.marginal.qtl(X.std, R) %>%
    rename(rs = snp, ensg = gene)  %>%
        mutate(rs = as.character(rs), ensg = as.character(ensg)) %>%
            left_join(genes, by = 'ensg') %>%
                left_join(x.bim %>% select(-chr), by = 'rs') %>%
                    filter((tss - cis.dist) < snp.loc & snp.loc < (tes + cis.dist)) %>%
                        filter(!ensg %in% rm.ensg)

write.mat(marginal.qtl, file = gzfile(out.qtl.file))
write.mat(samples.matched, file = gzfile(out.sample.file))
write.mat(genes, file = gzfile(out.gene.file))
write.mat(R, file = gzfile(out.resid.file))
