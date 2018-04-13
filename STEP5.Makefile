## Run FQTL

CHR := $(shell seq 1 22)

all: $(foreach chr, $(CHR), jobs/step5/run_$(chr)_fqtl.txt.gz jobs/step5/run_$(chr)_null.txt.gz)

jobs/step5/run_%_fqtl.txt.gz: 
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 processed/combined/chr$*/*.txt.gz | sed 's/.txt.gz//' | awk -vOFS=' ' -F'/' '{ print "./make.fqtl.R" OFS ($$0 ".txt.gz") OFS ("geno/" $$3) OFS ("result/fqtl/" $$3 "/" $$4) OFS "FALSE" }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N fqtl_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step5/run_%_null.txt.gz: 
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@ls -1 processed/combined/chr$*/*.txt.gz | sed 's/.txt.gz//' | awk -vOFS=' ' -F'/' '{ print "./make.fqtl.R" OFS ($$0 ".txt.gz") OFS ("geno/" $$3) OFS ("result/fqtl/" $$3 "/" $$4) OFS "TRUE" }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N null_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

