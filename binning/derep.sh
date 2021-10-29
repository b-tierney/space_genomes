#!/bin/bash
​
input_folder=$1

# merge bins into single folder

mkdir "${input_folder}"/all_bins

find "${input_folder}"/*_bins.checkm/bins/*/*fa > "${input_folder}"_binlocs




cp -r  "${input_folder}"/all_bins
​
# run drep to cluster assembled genomes from the same species
dRep dereplicate drep_workdir \
	-g "${input_folder}"/*.fa \
	-p 30 \
	-comp 50 \
	-pa 0.9 \
	-sa 0.95 \
	-nc 0.30 \
	-cm larger \
	--genomeInfo checkm_output.csv