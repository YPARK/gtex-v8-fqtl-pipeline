#!/bin/awk -f
BEGIN { 
    FS = "\t"
    print "chromosome" FS "position" FS "effect_allele" FS "non_effect_allele" FS "zscore" FS "pvalue"
}
$3 == CHR && $5 ~ /[ACTG]/ && $6 ~ /[ACTG]/ {
    printf $3 FS $4 FS $5 FS $6 FS;
    printf("%.4f" FS "%.4e\n", $10, $11)
}
