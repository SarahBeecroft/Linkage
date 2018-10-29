#!/bin/bash

DATA_DIR=/home/san/sbeecroft/old_linkdatagen

cd $DATA_DIR

for file in D18_0389 D18_0394 D18_0396
do 

nice /home/san/sbeecroft/atesta/home/san/atesta/src/samtools-0.1.19/samtools-0.1.19/samtools mpileup -d10000 -q13 -Q13 -gf /REFERENCE/Sequence/ucsc.hg19.fasta -l /home/san/sbeecroft/old_linkdatagen/annotHapMap2L.txt $DATA_DIR/$file.bam | nice /home/san/sbeecroft/atesta/home/san/atesta/src/samtools-0.1.19/bcftools/bcftools view -cg -t0.5 - > $file.HM.vcf

done
