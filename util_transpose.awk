BEGIN {
}
{ 
    for (i=1; i<=NF; i++)  {
	v=$i;
	gsub(/ /, "", v);
	v = length(v) > 0 ? $i : "NA";
        a[NR,i] = v;
    }
    p = NF > p ? NF : p;
    n = NR
}
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=n; i++){
            str=str FS a[i,j];
        }
        print str;
    }
}
