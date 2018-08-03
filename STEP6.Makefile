
BLK := 10

CHR := $(shell seq 1 22)

GWAS_DIR := gtex-gwas-hg38-imputed

GWAS_FILES := $(shell ls -1 $(GWAS_DIR)/*.txt.gz)

GWAS := $(shell cat $(GWAS_DIR)/gtex-gwas-share.tsv | tail -n+2 | cut -f2 | sed 's/.txt.gz//g')

all:

data: jobs/20180729/data.txt.gz

combine: jobs/20180729/combine.txt.gz

twas: $(foreach gwas, $(GWAS), jobs/20180729/$(gwas).twas.txt.gz)

twas-long: $(foreach gwas, $(GWAS), jobs/20180729/$(gwas).long-twas.txt.gz)

################################################################
jobs/20180729/data.txt.gz: $(foreach gf, $(GWAS_FILES), jobs/20180729/data-$(gf)-jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=2:00:00 -b y -j y -N GWAS_DATA -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
# % = gtex-gwas-hg38/UKB_20002_1261_self_reported_multiple_sclerosis.txt.gz
jobs/20180729/data-%-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 22 | awk '{ print "./make.distribute.gwas.sh" FS "$*" FS $$1 }' > $@

################################################################
# % = UKB_20002_1261_self_reported_multiple_sclerosis
jobs/20180729/%.twas.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 22 | awk -vB=$(BLK) -vGWAS=$* '{ CHR = $$1; (("zcat stat/fqtl_" CHR ".txt.gz | tail -n+2 | wc -l") | getline N); nb = int(N/B); for(b = 0; b <= nb; ++b) print "./make.twas.R" FS ("stat/fqtl_" CHR ".txt.gz") FS ("gwas_data/" GWAS "/" CHR ".txt.gz") FS ("geno/chr" CHR) FS ((B * b + 1) ":" (B* (1 + b))) FS ("twas/" GWAS "/chr" CHR "_" b ".txt.gz") }' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N TWAS_$(shell echo $* | sed 's/\//_/g') -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180729/%.long-twas.txt.gz: jobs/20180729/%.twas.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N LONG_$(shell echo $* | sed 's/\//_/g') -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
jobs/20180729/combine.txt.gz: $(foreach gwas, $(GWAS), jobs/20180729/combine-$(gwas).txt)
	cat $^ | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N TWAS_COMBINE -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180729/combine-%.txt: 
	echo "./make.combine_twas.R twas/$* twas/$*.txt.gz" > $@



