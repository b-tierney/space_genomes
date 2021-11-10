#!/bin/bash

# run snippy core

files=$(cat "space_genomes_mapped" | sed "s/$/_"$2"/g" | tr '\n' '\t')

snippy-core --ref "$1" --prefix="$2" $files