#!/bin/bash

echo '''
Used for making the HM.vcf files for the linkdatagen vcf2linkdatagen.pl and linkdatagen.pl programs. Requires space seperated list of sample names.

Usage:
./gen_HM.vcf.sh <data_dir> <samtools_dir> <reference> <sample_list> <annot_dir>

Cheat: if your current working directory is what you would need to put in one of the above arguments, you can use . instead of typing out the full name e.g. 
./gen_HM.vcf.sh . . hg19.fa sample_list.txt .

INFO:

<data_dir> is the full path where your bam files are stored (must all be in same directory)
<samtools_dir> is the full path of the samtools-0.1.19 directory (e.g. /data/linkage/samtools-0.1.19)
<reference> is where your reference genome fasta file is located, plus the name of the reference file (i.e. /data/linkage/ucsc.hg19.fasta)
<sample_list> is your list of sample IDs (one sample on each line), which identify your bam files. E.g. D17-0323.bam and D98-0089.bam would be listed as simply D17-0323 D98-0089
<annot_dir> is where your linkdatagen reference files are located (e.g. /data/linkage)
'''

DATA_DIR=$1 
samtools_dir=$2
reference=$3
sample_list=$4
annot_dir=$5

echo 'DATA_DIR='$DATA_DIR
echo 'samtools_dir='$samtools_dir
echo 'reference='$reference
echo 'sample_list='$sample_list
echo 'annot_dir='$annot_dir

cd $DATA_DIR

while IFS= read -r file
do
echo 'processing ' $file
$samtools_dir/samtools-0.1.19/samtools mpileup -d10000 -q13 -Q13 -gf $reference -l $annot_dir/annotHapMap2U.txt $DATA_DIR/$file.bam | $samtools_dir/samtools-0.1.19/bcftools/bcftools view -cg -t0.5 - > $DATA_DIR/$file.HM.vcf
done < "$sample_list"
