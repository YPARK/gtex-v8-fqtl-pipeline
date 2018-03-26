
V8READ := rawdata/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count.gct.gz
GENCODE	:= rawdata/references/gencode.v26.GRCh38.genes.gtf

################################################################
## 2. Break down large expression matrix into small pieces:
##    chr x tissue = 22 x 49

CHR := $(shell seq 1 22 | awk '{ print "chr" $$1 }')

all : $(foreach tis, $(shell cat data/tissues.txt | cut -f1), \
                     data/$(tis)/samples.txt.gz \
                     data/$(tis)/rsem-count.txt.gz \
                     data/$(tis)/vst.txt.gz data/$(tis)/sz.txt.gz)

# % = $(tis)
data/%/samples.txt.gz: data/covariates.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat data/covariates.txt.gz | awk -vTIS=$* -F'\t' '$$(NF) == TIS' | gzip > $@

# % = $(tis)
data/%/rsem-count.txt.gz: data/coding.genes.txt.gz data/%/samples.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $(V8READ) | awk -F'\t' -vROWS=$$(zcat $< | awk '{ printf ((n++ > 0)? "," : "") $$NF }') -f util_subset_rows.awk | awk -F'\t' -vCOLS=$$(zcat data/$*/samples.txt.gz | awk '{ print $$3 }' | awk '{ printf ((n++ > 0)? "," : "") $$NF }') -f util_subset_cols.awk | awk -F'\t' -f util_transpose.awk | gzip > $@

# % = $(tis)
data/%/vst.txt.gz: data/%/rsem-count.txt.gz data/%/samples.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh ./make.preprocess.R data/$*

data/%/sz.txt.gz: data/%/vst.txt.gz
