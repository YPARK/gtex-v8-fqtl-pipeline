#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

library(dplyr)
library(tidyr)
options(stringsAsFactors = FALSE)

rnaseq.sample.file <- argv[1] # e.g., rnaseq.sample.file = 'data/rnaseq.samples.txt'
samp.annot.file <- argv[2]  # e.g., samp.annot.file = 'rawdata/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'
subj.annot.file <- argv[3] # e.g., subj.annot.file = 'rawdata/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt'
out.file <- argv[4] # e.g., out.file = 'data/covariates.txt.gz'

.read.table <- function(...) read.table(..., quote = '', comment.char = '', sep = '\t',
                                        stringsAsFactors = FALSE, fill = TRUE)

strip.br <- function(...) gsub(pattern = '[[:punct:]]', replacement = '', ...)

strip.ws <- function(...) gsub(pattern = '\\s+', replacement = '_', ...)

clean.str <- function(s) strip.ws(strip.br(s))

.scale <- function(...) as.numeric(scale(...))

samples <- .read.table(rnaseq.sample.file)
colnames(samples) <- c('SAMPID', 'SUBJID', 'data.pos')

samp.tab <- .read.table(samp.annot.file, header = TRUE) %>%
    select(SAMPID, SMRIN, SMTSISCH, SMTSD) %>%
        mutate(tis.name = sapply(SMTSD, clean.str))


subj.tab <- .read.table(subj.annot.file, header = TRUE) %>%
    mutate(SUBJID = sapply(SUBJID, function(s) strsplit(s, '[-]')[[1]][2])) %>%
        select(SUBJID, SEX, AGE, HGHT, WGHT, BMI) %>%
            mutate(SEX.std = .scale(SEX),
                   AGE.std = .scale(AGE),
                   HGHT.std = .scale(HGHT),
                   WGHT.std = .scale(WGHT),
                   BMI.std = .scale(BMI))

out.tab <- samples %>% left_join(samp.tab, by = 'SAMPID') %>%
    left_join(subj.tab, by = 'SUBJID')

valid.tis.tab <- out.tab %>% group_by(tis.name) %>%
    summarize(ntis = length(unique(SAMPID))) %>%
        filter(ntis >= 100) %>%
            arrange(tis.name) %>% mutate(tis.idx = 1:n()) %>%
                mutate(tis.dir = paste(tis.idx, tis.name, sep = '_'))

out.tab <- out.tab %>% filter(tis.name %in% valid.tis.tab$tis.name) %>%
    select(-SMTSD) %>% left_join(valid.tis.tab, by = 'tis.name')

write.table(out.tab, file = gzfile(out.file), sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = TRUE)
