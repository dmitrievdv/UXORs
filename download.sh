#!/bin/bash
input="UXORs.data" 
while read -r line; do
    words=( $line )
    name=${words[0]}
    RA=${words[1]}
    DEC=${words[2]}
    URL="https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=$RA&dec=$DEC&y=15&x=15&units=px&product=SPOC&sector=All"
    if [ ! -d "$name" ]
    then
        mkdir $name
        curl -0 "$URL" -o "$name/$name-cut.zip"
        cd $name
        unzip $name-cut.zip
        rm -rf $name-cut.zip
        cd ..
    fi
done < $input