#!/bin/bash

#could also use just snippy-multi to the same effect

snippy --cpus 5 --outdir "$1"_"$5" --ref $4 --R1 $2 --R2 $3