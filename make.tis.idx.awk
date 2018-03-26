#!/bin/gawk -f

BEGIN {
    FS = "\t"
}

{
    if(!($2 in tis2idx)) tis2idx[$2] = ++idx
    rows[NR] = $0
    row2tis[NR] = $2
}

END {
    for(r in rows) {
	t = row2tis[r]
	print rows[r] FS tis2idx[t]
    }
}
