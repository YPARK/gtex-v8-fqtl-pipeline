## confoundeer correction for cis-eqtl calling

all: jobs/gene_segments.txt \
    $(foreach tis, $(shell cat data/tissues.txt | cut -f1), \
     jobs/step3/$(tis)-confounder.txt.gz)

long: $(foreach tis, $(shell cat data/tissues.txt | cut -f1), jobs/step3/$(tis)-confounder-long.txt.gz)

jobs/gene_segments.txt: make.job.segments.R data/coding.genes.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@Rscript --vanilla $^ $@

jobs/step3/%-confounder.txt.gz: jobs/gene_segments.txt
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ -d processed/expression/$* ] || mkdir -p processed/expression/$*
	@cat $< | awk '{ print "./make.cis-correction.R" FS "$*" FS $$1 FS $$2 ":" $$3 FS "processed/expression/$*/nb_" $$2 ".txt.gz"; }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N conf_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3/%-confounder-long.txt.gz: jobs/step3/%-confounder.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ -d processed/expression/$* ] || mkdir -p processed/expression/$*
	@zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=16:00:00 -b y -j y -N long_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@


################################################################
## Utilities
PLINKZIP := https://www.cog-genomics.org/static/bin/plink170906/plink_linux_x86_64.zip

bin/plink:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	curl $(PLINKZIP) -o bin/plink.zip
	unzip bin/plink.zip -d bin/

