#!/bin/bash -l

if [ $# -lt 2 ]; then
    echo "Need \$1 = gwas_file"
    echo "Need \$2 = chromosome_number"
    exit 1
fi

gwas_file=$1 # e.g., gtex-gwas-hg38/IBD.EUR.Crohns_Disease.txt.gz
chr=$2       # e.g., chr=21

# simply break down gwas chromosome by chromosome
out_dir=gwas_data/$(basename $gwas_file .txt.gz)/
mkdir -p $out_dir
out_file=$out_dir/$chr.txt.gz
[ -f $out_file ] || \
    zcat $gwas_file | \
	./make.distribute.gwas.awk -vCHR=chr$chr | \
	gzip > $out_file

exit 0
