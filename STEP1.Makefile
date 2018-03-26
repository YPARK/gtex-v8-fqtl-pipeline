
V8READ := rawdata/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_expected_count.gct.gz
SAMP := rawdata/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt
SUBJ := rawdata/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt
GENCODE	:= rawdata/references/gencode.v26.GRCh38.genes.gtf

META := data/coding.genes.txt.gz \
  data/covariates.txt.gz \
  data/rnaseq.samples.txt \
  data/tissues.txt

## 1. Take coding genes and clean up covariate matrices

all : $(META)

data/rnaseq.samples.txt:
	zcat $(V8READ) | head -n3 | tail -n1 | tr '\t' '\n' | awk -F'\t' '{ split($$1,name,"-"); print $$1 FS name[2] FS NR }' | tail -n+3 > $@

data/coding.genes.txt.gz: ./make.coding.genes.R $(GENCODE) data/rnaseq.rows.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

data/rnaseq.rows.txt.gz: $(V8READ)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | cut -f 1 | awk -F'\t' '{ print $$1 FS NR }' | gzip > $@

data/covariates.txt.gz: make.samples.cov.R data/rnaseq.samples.txt $(SAMP) $(SUBJ)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

data/tissues.txt: data/covariates.txt.gz
	zcat $< | tail -n+2 | awk -F'\t' '{ out[$$18] = ($$19 FS $$6 FS $$17 FS $$6) } END { for(j in out) print j FS out[j] }' | sort -k1n | cut -f 2- > $@

clean:
	rm -f $(META)
