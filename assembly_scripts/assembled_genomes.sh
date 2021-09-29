#!/bin/bash

basename=$(echo $1 | cut -f1-3 -d_)

spades.py -1 $1 -2 $2 -o "$basename"_assembled -t 20 -m 200

reformat.sh in="$basename"_assembled/scaffolds.fasta out="$basename"_assembled/scaffolds_filtered.fasta minlength=1000

quast.py -o "$basename"_assembled/ "$basename"_assembled/scaffolds_filtered.fasta