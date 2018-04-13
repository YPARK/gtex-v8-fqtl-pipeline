## simply combine confounder-corrected expressions gene by gene

all: jobs/step4/combine.txt.gz

jobs/step4/combine.txt.gz: jobs/gene_segments.txt
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | awk '{ print "./make.combine_tissue.R" FS $$2 FS "processed/combined/" $$1 }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=1g -l h_rt=1:00:00 -b y -j y -N combine_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@
