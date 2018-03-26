#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
library(readr)
library(DESeq2)
library(vsn)
library(dplyr)

data.dir <- argv[1] # e.g., data.dir = 'data/35_Ovary'

data.file <- data.dir %&&% '/rsem-count.txt.gz'
sample.file <- data.dir %&&% '/samples.txt.gz'
vst.file <- data.dir %&&% '/vst.txt.gz'
sz.file <- data.dir %&&% '/sz.txt.gz'

remove.zero <- function(Y) {
    Y[ Y < 1 ] <- NA
    return(Y)    
}

read.tsv <- function(...) {
    readr::read_tsv(..., col_names = FALSE) %>%
        as.matrix()
}

adjust.size <- function(mat, sf) {
    sweep(mat, 1, sf, `/`)
}

Y <- read.tsv(data.file)

plt.file.0 <- data.dir %&&% '/mean_sd_plot_Y.pdf'
pdf(file = plt.file.0, useDingbats = FALSE)
meanSdPlot(t(Y))
dev.off()

plt.file.0 <- data.dir %&&% '/mean_sd_plot_log2Y.pdf'
pdf(file = plt.file.0, useDingbats = FALSE)
meanSdPlot(log2(1/2 + t(Y)))
dev.off()

covar.cols <- c('SAMPID', 'SUBJID', 'data.pos', 'SMRIN', 'SMTSISCH',
                'tis.name', 'SEX', 'AGE', 'HGHT', 'WGHT', 'BMI', 'SEX.std', 'AGE.std',
                'HGHT.std', 'WGHT.std', 'BMI.std', 'ntis', 'tis.idx', 'tis.dir')

covar.tab <- readr::read_tsv(sample.file, col_names = covar.cols)

covar.tab.sub <- covar.tab %>% dplyr::select(dplyr::ends_with('.std')) %>%
    as.data.frame()

if(var(covar.tab.sub$SEX.std) < 1e-4) {
    dds <- DESeqDataSetFromMatrix(t(round(Y)), covar.tab.sub,
                                  ~ AGE.std + HGHT.std + WGHT.std + BMI.std + 1)
} else {
    dds <- DESeqDataSetFromMatrix(t(round(Y)), covar.tab.sub,
                                  ~ SEX.std + AGE.std + HGHT.std + WGHT.std + BMI.std + 1)
}

dds <- estimateSizeFactors(dds)

dds <- estimateDispersions(dds)

vst.out <- varianceStabilizingTransformation(dds, blind = FALSE)

Y.vst <- assay(vst.out) %>% t() %>% scale()

plt.file.1 <- data.dir %&&% '/mean_sd_plot_Y_vst.pdf'
pdf(file = plt.file.1, useDingbats = FALSE)
meanSdPlot(t(Y.vst))
dev.off()

plt.file.1 <- data.dir %&&% '/mean_sd_plot_log2Y_vst.pdf'
pdf(file = plt.file.1, useDingbats = FALSE)
meanSdPlot(log2(1/2 + t(Y.vst)))
dev.off()

write_tsv(data.frame(round(Y.vst, 4)), path = vst.file, col_names = FALSE)

write_tsv(data.frame(round(sizeFactors(dds), 4)), path = sz.file, col_names = FALSE)

