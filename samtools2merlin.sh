#!/bin/bash

###Usage: Edit the DATA_DIR variable to the directory where your bams,  other input files, linkdatagen references, and linkage scripts are stored. Edit the sample/file names below to your files of interest. Save and close the script.

#'file' variable refers to the sample ID you've given to your input files. e.g. OPMD. input files should be $file.vcflist.txt, $file.ws, $file.brlmm, $file.ped

DATA_DIR=/home/san/sbeecroft/old_linkdatagen
cd $DATA_DIR
echo -en "Enter a space seperated list of your sample IDs and hit return: \n"
read sample
echo -en $sample" have been accepted as sample IDs. \n"
echo -en "Enter input file prefix and hit return: \n"
read file
echo -en "File prefix is "$file" \n"
echo -en "Enter inheritance model (options are dominant or recessive- make sure to check your spelling!): \n "
read model
echo -en $model" has been entered as inheritance pattern. \n"

for sample in $(echo $sample)
do 

nice /home/san/sbeecroft/atesta/home/san/atesta/src/samtools-0.1.19/samtools-0.1.19/samtools mpileup -d10000 -q13 -Q13 -gf /REFERENCE/Sequence/ucsc.hg19.fasta -l $DATA_DIR/annotHapMap2L.txt $DATA_DIR/$sample.bam | nice /home/san/sbeecroft/atesta/home/san/atesta/src/samtools-0.1.19/bcftools/bcftools view -cg -t0.5 - > $DATA_DIR/$sample.HM.vcf

done

cd $DATA_DIR
nice perl vcf2linkdatagen.pl -annotfile annotHapMap2.txt -pop CEU -mindepth 20 -missingness 0 -idlist $file.vcflist.txt > $file.brlmm
nice perl linkdatagen.pl -data m -pedfile $file.ped -whichSamplesFile $file.ws -callFile $file.brlmm -annotFile annotHapMap2.txt -pop CEU -binsize 0.3 -mindist 0.15 -prog me -outputDir $file.HapMap2
cd /home/san/sbeecroft/old_linkdatagen/$file.HapMap2_merlin/genome/
cp /home/san/sbeecroft/merlin-1.1.2/executables/$model.param.tbl temp.tbl
nice merlin -d merlin_autosome_$file.dat -p merlin_genome_$file.ped -f merlin_genome_$file.freq -m merlin_genome_$file.map --prefix merlin_autosome_$file --error
nice pedwipe -d merlin_autosome_$file.dat -p merlin_genome_$file.ped -e merlin_autosome_$file.err
mv wiped.dat merlin_autosome_wiped_$file.dat
mv wiped.ped merlin_autosome_wiped_$file.ped
nice minx -d merlin_X_$file.dat -p merlin_genome_$file.ped -f merlin_genome_$file.freq -m merlin_genome_$file.map --prefix merlin_X_$file --error
nice pedwipe -d merlin_X_$file.dat -p merlin_genome_$file.ped -e merlin_X_$file.err
mv wiped.dat merlin_X_wiped_$file.dat
mv wiped.ped merlin_X_wiped_$file.ped
rm wiped.freq
mv temp.tbl $model.param.tbl
nice merlin -d merlin_autosome_wiped_$file.dat -p merlin_autosome_wiped_$file.ped -f merlin_genome_$file.freq -m merlin_genome_$file.map --smallswap --megabytes:9999 --founders --MarkerNames --pairs --ibd --extended --best --exp --model $model.param.tbl --pdf --tabulate --prefix merlin_autosome_$file
nice minx -d merlin_X_wiped_$file.dat -p merlin_X_wiped_$file.ped -f merlin_genome_$file.freq -m merlin_genome_$file.map --smallswap --megabytes:9999 --founders --MarkerNames --pairs --ibd --extended --best --exp --model $model.param.tbl --pdf --tabulate --prefix merlin_X_$file
