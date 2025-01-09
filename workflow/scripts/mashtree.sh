#!/bin/bash

MASHTREE_INPUT="$1"
MASHTREE_TEMP_DIR="$2"
mkdir -p "$MASHTREE_TEMP_DIR"
while read line
    do
        if [[ ${line:0:1} == '>' ]]
        then
            outfile="$MASHTREE_TEMP_DIR"/${line#>}.fa
            echo $line > $outfile
        else
            echo $line >> $outfile
        fi
    done < "$MASHTREE_INPUT"
