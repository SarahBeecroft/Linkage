#! /bin/bash

for file in maori_homozygosity_mapping

do

cd /home/san/sbeecroft/old_linkdatagen/

nice perl vcf2linkdatagen.pl -annotfile annotHapMap2.txt -pop CEU -mindepth 20 -missingness 0 -idlist $file.vcflist.txt > $file.brlmm \
&& nice perl linkdatagen.pl -data m -pedfile $file.ped -whichSamplesFile $file.ws -callFile $file.brlmm -annotFile annotHapMap2.txt -pop CEU -binsize 0.3 -mindist 0.15 -prog cp -outputDir $file.HapMap2 \

done
