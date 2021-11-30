#!/bin/bash

pyseer --lmm --phenotypes pyseer_metadata.tsv --pres $1 --similarity distances_with_bins > pyseer_"$1"
