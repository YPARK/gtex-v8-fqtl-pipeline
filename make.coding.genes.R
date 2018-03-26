#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) {
    q()
}

options(stringsAsFactors = FALSE)
library(dplyr)
library(rtracklayer)

## select coding genes in GTF
gtf.file <- argv[1]       # e.g., gtf.file <- 'rawdata/references/gencode.v26.GRCh38.genes.gtf'
data.gene.file <- argv[2] # e.g., data.gene.file <- 'data/rnaseq.genes.txt.gz'
out.file <- argv[3]

gtf.tab <- readGFF(gtf.file, tags = c('gene_id', 'gene_name', 'transcript_name', 'gene_type'))

coding.genes <- gtf.tab %>% mutate(chr = seqid, ensg = gene_id) %>%
    filter(gene_type == 'protein_coding', type == 'transcript') %>%
        select(chr, start, end, strand, ensg, gene_name) %>%
            arrange(chr, start)

data.genes <- read.table(data.gene.file)
colnames(data.genes) <- c('ensg', 'data.loc')

chr.names <- paste('chr',c(1:22),sep='')

out.tab <- coding.genes %>%
    left_join(data.genes, by = 'ensg') %>%
        na.omit() %>%
            filter(chr %in% chr.names)

write.table(out.tab, file = gzfile(out.file), sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
