#!/bin/gawk -f

BEGIN {
    FS = "\t";
    OFS = "\t";
}

{
    s = $1;
    tis = $14;

    gsub(/[()]+/,"", tis);
    gsub(/[ -]+/, "_", tis);

    samp2tis[s] = tis;
    ntis[tis] ++;

}

END {

    for(i in samp2tis) {
	t = samp2tis[i]
	if(ntis[t] >= 100)
	    print i FS t FS ntis[t]
    }
}
