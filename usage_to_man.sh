#!/bin/bash

#Converts the Fastaq python scripts usage into man pages
#The man pages are placed in the man folder of the main Fastaq directory  

source='scripts'
destination='man'

for i in $source/fastaq_*; do
    j=$(echo $i | sed "s/$source\///");
    out_file=$destination/$j.1
    help2man -m $j -n $j --no-discard-stderr $i |
    grep -v $j': error: too few arguments' |
    sed "s/usage://gi" > $out_file;
done



for file in $destination/*.1; do

    out_file=$(echo $file | sed 's/\.1/\.2/');
    filename=$(echo $file | sed "s/\($destination\/*\)\([a-zA-Z_0-9]* *\)\.1/\2/")
    while read line
    do
	echo $line |
	sed "s/\-\s\($filename *\)/\1/" |
	sed "s/\(.TH *\)\([a-zA-Z_0-9]* *\)\(\"[0-9\"]* *\)\(\"[a-zA-Z0-9_ ]*\" *\)\(\"[a-zA-Z0-9_\<\>\/ ]*\" *\)\(\"[a-zA-Z0-9_]*\" *\)/\1 \2 \6 \3/" |
	sed "s/\(The\sfull\sdocumentation\sfor*\)/$filename -h will display the original python documentation/" |
	sed "s/^\.B//" |
	sed "s/^is\smaintained\sas\sa\sTexinfo\smanual.\sIf\sthe//" |
	sed "s/^\.B//" |
	sed "s/\sinfo//" |
	sed "s/^and//" |
	sed "s/^programs\sare\sproperly\sinstalled\sat\syour\ssite\,\sthe\scommand//" |
	sed "s/\(.IP\r*\)//" |
	sed "s/^should\sgive\syou\saccess\sto\sthe\scomplete\smanual\.//"

    done < $file > $out_file
cat <<EOF >> $out_file
.SH "AUTHOR"
.sp
$filename was originally written by Martin Hunt (mh12@sanger\&.ac\&.uk)
.SH "COPYING"
.sp
Wellcome Trust Sanger Institute Copyright \(co 2013 Wellcome Trust Sanger Institute This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version\&.
EOF

    mv $out_file $file


done
