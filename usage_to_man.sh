#!/bin/bash
for i in scripts/fastaq_*; do
    j=$(echo $i | sed 's/scripts\///');
    out_file=temp_man/$j.1
    help2man -m $j -n $j --no-discard-stderr $i |
    grep -v $j': error: too few arguments' |
    sed "s/usage://gi" > $out_file;
done


for file in temp_man/*.1; do

    out_file=$(echo $file | sed 's/\.1/\.2/');
    while read line
    do

	echo $line | sed "s/\(.TH *\)\([a-zA-Z_0-9]* *\)\(\"[0-9\"]* *\)\(\"[a-zA-Z0-9_ ]*\" *\)\(\"[a-zA-Z0-9_\<\>\/ ]*\" *\)\(\"[a-zA-Z0-9_]*\" *\)/\1 \2 \3 \4/"

    done < $file > $out_file
    echo $file
    echo $out_file
    `mv "$outfile" "$file"`

done
