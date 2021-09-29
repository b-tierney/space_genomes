#!/bin/bash

mash sketch -s 10000 -o mash_sketch bacterial_genomes/*.fa
mash dist mash_sketch.msh mash_sketch.msh > apit_mash.tsv
