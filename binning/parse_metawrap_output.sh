#!/bin/bash

# parse output of metawrap

#binning output made by ls */*assembled/scaffolds.fasta_binning/ | grep : | cut -f1 -d: > binning_output
while read p; do
	mkdir -p "${p}"/forderep
	for file in $(ls "${p}"/*bins/*fa); do
			bintype=$(echo $file | rev | cut -f2 -d/ | rev);
			filename=$(echo $file | rev | cut -f1 -d/ | rev);
			echo cp "${p}" "${p}"/forderep/"${bintype}_${filename}";
	done
done<binning_output