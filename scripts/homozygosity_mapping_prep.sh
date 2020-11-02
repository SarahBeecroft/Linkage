#!/bin/bash

echo 'usage is ./homozygosity_mapping_prep.sh <file_prefix> (e.g. OPDM)'

file=$1

perl vcf2linkdatagen.pl -annotfile annotHapMap2U.txt -pop CEU -mindepth 20 -missingness 0 -idlist $file.vcflist.txt > $file.brlmm

perl linkdatagen.pl -data m -pedfile $file.ped -whichSamplesFile $file.ws -callFile $file.brlmm -annotFile annotHapMap2U.txt -pop CEU -binsize 0.3 -mindist 0.15 -prog cp -outputDir $file.HapMap2
