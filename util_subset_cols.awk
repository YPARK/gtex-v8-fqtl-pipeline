#!/bin/awk -f

BEGIN{
    split(COLS, cols, ",");
    for(j=1; j<=length(cols); ++j) {
	col2pos[j] = cols[j];
	pos2col[cols[j]] = j;
    }
    IGNORECASE = 1;
}
{
    c = col2pos[1];
    v = $c;
    if($c ~ /na/ || length($c) == 0){
	v = "NA"
    }
    printf v

    for(j=2; j<=length(cols); ++j){
	c = col2pos[j];
	v = $c;
	if($c ~ /na/){
	    v = "NA"
	}
	printf FS v;
    }
    printf "\n";
}

