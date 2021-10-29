#!/bin/bash

conda activate p2

metawrap bin_refinement -o "${1}"_refined_bins -t 5 -A "${1}"/metabat2_bins -B /concoct_bins -C "${1}"/maxbin2_bins