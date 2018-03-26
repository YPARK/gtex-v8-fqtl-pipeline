#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) {
    q()
}

## What this script does:
## 1. Find non-genetic confounder variables out of the chromosome
## 2. Generate confounder-corrected expression
## 3. Remove potential confounders

tis.name <- argv[1]     # e.g., tis.name = '20_Breast_Mammary_Tissue'
chr <- argv[2]          # e.g., chr = 'chr1'
gene.idx.str <- argv[3] # e.g., gene.idx.str = '1:20'
out.file <- argv[4]     # e.g., out.file = 'temp.gz'

options(stringsAsFactors = FALSE)
source('Util.R')

if(file.exists(out.file)) {
    log.msg('File %s already exists\n', out.file)
    q()
}


library(fqtl)
library(dplyr)
library(readr)
library(methods)

dir.create(dirname(out.file), recursive = TRUE)

#########################
## important parameter ##
#########################

## We will not consider confouders that are correlated with
## cis-regulatory genotypes

ctrl.degree <- 10
p.val.cutoff <- 1e-4

################################################################
gene.idx <- eval(parse(text = gene.idx.str))
gene.file <- 'data/coding.genes.txt.gz'

temp.dir <- system('mkdir -p /broad/hptmp/ypp/gtexv8/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/gtexv8/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

gene.cols <- c('chr', 'lb', 'ub', 'strand', 'ensg', 'hgnc', 'remove')
covar.cols <- c('SAMPID', 'SUBJID', 'data.pos', 'SMRIN', 'SMTSISCH',
                'tis.name', 'SEX', 'AGE', 'HGHT', 'WGHT', 'BMI', 'SEX.std', 'AGE.std',
                'HGHT.std', 'WGHT.std', 'BMI.std', 'ntis', 'tis.idx', 'tis.dir')

data.file <- 'data/' %&&% tis.name %&&% '/rsem-count.txt.gz'
sample.file <- 'data/' %&&% tis.name %&&% '/samples.txt.gz'
sz.file <- 'data/' %&&% tis.name %&&% '/sz.txt.gz'

genes.info <- read_tsv(gene.file, col_names = gene.cols) %>%
    dplyr::select(-remove)

covar.tab <- readr::read_tsv(sample.file, col_names = covar.cols)

size.factors <- read_tsv(sz.file, col_names = 'sz')

################################################################
Y <- read_tsv(data.file, col_names = FALSE) %>%
    adjust.size.factor() %>%
        stdize.count() %>%
            as.matrix()

Y1 <- Y %c% gene.idx
genes.Y1 <- genes.info %r% gene.idx

genes.bound.y1 <- genes.Y1 %>%
    group_by(chr) %>%
        summarize(lb.min = min(lb), ub.max = max(ub))

## out-of-chromosome or out-of-regulatory
other.genes <- genes.info %>% mutate(idx = 1:n()) %>%
    dplyr::filter(!(chr %in% genes.Y1$chr))

other.genes.cis <-
    genes.info %>% mutate(idx = 1:n()) %>%
        dplyr::filter(chr %in% genes.Y1$chr) %>%
            left_join(genes.bound.y1) %>%
                dplyr::filter(ub < lb.min | lb > ub.max) %>%
                    select(-lb.min, -ub.max)

other.genes <- bind_rows(other.genes, other.genes.cis)

other.idx <- other.genes$idx
Y0 <- Y %c% other.idx

## remove genes with too many zeros or NAs
y1.rm <- gene.filter.nb(Y1)
y0.rm <- gene.filter.nb(Y0)

if(length(y1.rm) > 0) {
    Y1 <- Y1 %c% -y1.rm
    genes.Y1 <- genes.Y1 %r% -y1.rm
    gene.idx <- gene.idx[-y1.rm]
}

if(length(y0.rm) > 0) {
    Y0 <- Y0 %c% -y0.rm
}

if(nrow(genes.Y1) < 1) {
    write_tsv(x = data.frame(), path = out.file)
    q()
}

################################################################
nn.ctrl <- min(ctrl.degree, ncol(Y0))

log.msg('Find %d control genes for each target gene\n', nn.ctrl)

ctrl.idx <- find.cor.idx(Y1 %>% trans.normal(),
                         Y0 %>% trans.normal(),
                         n.ctrl = nn.ctrl)

Y0.ctrl <- Y0 %c% ctrl.idx %>% as.matrix()
genes.ctrl <- other.genes %r% ctrl.idx

################################################################
## take cis-regulatory genotypes of Y1
cis.dist <- 1e6
lb <- pmax(min(genes.Y1$lb) - cis.dist, 0)
ub <- pmax(max(genes.Y1$ub) + cis.dist, 0)

plink <- subset.plink(glue('geno/', chr), chr, lb, ub, temp.dir)

system('rm -r ' %&&% temp.dir)

x.fam <- plink$FAM %>%
    dplyr::select(iid) %>%
        dplyr::rename(SUBJID = iid) %>%
            mutate(x.pos = 1:n())

x.bim <- plink$BIM

covar.matched <- covar.tab %>% mutate(y.pos = 1:n()) %>%
    left_join(x.fam) %>% na.omit()

x.pos <- covar.matched$x.pos
y.pos <- covar.matched$y.pos

################################################################
## just remove Y0 genes correlated with genotype of Y1
xx.std <- plink$BED %r% x.pos %>% scale() %>% rm.na.zero()

Y0.std <- Y0.ctrl %>% trans.normal() %r% y.pos %>% scale() %>% rm.na.zero()
colnames(Y0.std) <- 1:ncol(Y0.ctrl)

x.y0.qtl <- calc.qtl.stat(xx.std, Y0.std)

max.qtl <- x.y0.qtl %>%
    group_by(y.col)%>%
        slice(which.min(p.val))

rm.cols <- max.qtl %>%
    dplyr::filter(p.val < p.val.cutoff) %>%
        dplyr::select(y.col) %>%
            unique() %>%
                unlist(use.names = FALSE)

################################################################
opt.reg <- list(vbiter = 3000, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi.ub = -1, pi.lb = -4, tau = -4, do.hyper = TRUE,
                jitter = 0.1, model = 'nb', out.residual = TRUE,
                print.interv = 100,
                k = min(ncol(Y0.ctrl), ncol(Y1)))

covar.mat <- covar.tab %>%
    dplyr::select(dplyr::ends_with('.std')) %>%
        as.matrix() %>%
            scale() %>%
                as.data.frame() %>%
                    mutate(intercept = 1) %>%
                        as.matrix()

if(length(rm.cols) > 0){
    genes.ctrl <- genes.ctrl %r% (-rm.cols)
    Y0.ctrl <- Y0.ctrl %c% (-rm.cols)
}

log.msg('Identified %d Y0 covariates; removed %d genes\n',
        ncol(Y0.ctrl),
        length(rm.cols))

################################################################
if(ncol(Y0.ctrl) > 0) {

    y0.covar <- Y0.ctrl %>% trans.normal() %>% scale() %>% rm.na.zero()

    ## remove non-genetic confounders
    y1.out <- fqtl.regress(y = Y1, 
                           x.mean = y0.covar,
                           c.mean = covar.mat,
                           factored = TRUE,
                           options = opt.reg)

} else {

    ## remove non-genetic confounders
    y1.out <- fqtl.regress(y = Y1,
                           x.mean = covar.mat,
                           options = opt.reg)

}

ret <- y1.out$resid$theta %>% scale() %>% signif(digits = 4) %>% as.data.frame()
colnames(ret) <- gene.idx

write_tsv(x = ret, path = gzfile(out.file), col_names = TRUE)

log.msg('Successfully finished!\n')

