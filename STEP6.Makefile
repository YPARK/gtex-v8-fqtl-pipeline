
BLK := 200

CHR := $(shell seq 1 22)

GWAS_DIR := gtex-gwas-hg38

GWAS := $(shell cat $(GWAS_DIR)/gtex-gwas-share.tsv | tail -n+2 | cut -f2 | sed 's/.txt.gz//g')

all:

queue: $(foreach gwas, $(GWAS), jobs/step6/$(gwas).twas.txt.gz)

long: $(foreach gwas, $(GWAS), jobs/step6/$(gwas).long-twas.txt.gz)

# % = UKB_20002_1261_self_reported_multiple_sclerosis
jobs/step6/%.twas.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	for((chr=1;chr<=22;++chr)); do awk -vB=$(BLK) -vCHR=$${chr} -vN=$$(zcat stat/fqtl_$${chr}.txt.gz | tail -n+2 | wc -l) -vGWAS=$* 'BEGIN { nb = int(N/B); for(b = 0; b <= nb; ++b) print "./make.twas.R" FS ("stat/fqtl_" CHR ".txt.gz") FS ("$(GWAS_DIR)/" GWAS ".txt.gz") FS ("geno/chr" CHR) FS ((B * b + 1) ":" (B* (1 + b))) FS ("twas/" GWAS "/chr" CHR "_" b ".txt.gz") }' | gzip >> $@; done
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=1:00:00 -b y -j y -N TWAS_$(shell echo $* | sed 's/\//_/g') -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step6/%.long-twas.txt.gz: jobs/step6/%.twas.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=12:00:00 -b y -j y -N LONG_$(shell echo $* | sed 's/\//_/g') -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

