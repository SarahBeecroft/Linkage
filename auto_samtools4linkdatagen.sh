#!/bin/bash

###Usage: Edit the DATA_DIR variable to the directory where your bams,  other input files, linkdatagen references, and linkage scripts are stored. Edit the sample/file names below to your files of interest. Save and close the script. Then, run using ./auto_samtools4linkdatagen.sh this may take some time to run.

DATA_DIR=/home/san/sbeecroft/old_linkdatagen
cd $DATA_DIR

for file in D18_0389 D18_0394 D18_0396
do 

nice /home/san/sbeecroft/atesta/home/san/atesta/src/samtools-0.1.19/samtools-0.1.19/samtools mpileup -d10000 -q13 -Q13 -gf /REFERENCE/Sequence/ucsc.hg19.fasta -l $DATA_DIR/annotHapMap2L.txt $DATA_DIR/$file.bam | nice /home/san/sbeecroft/atesta/home/san/atesta/src/samtools-0.1.19/bcftools/bcftools view -cg -t0.5 - > $DATA_DIR/$file.HM.vcf

done
