#!/bin/bash

echo "This is the correct script to use if you already have generated your HM.vcf files, ped, ws, and vcflist files. This script expects all inputs to be in the working directory, with Merlin binaries already installed to your PATH. If you need to generate the HM.vcfs, use the auto_samtools.sh script."
echo
echo "Enter the base name for your input files (i.e. family prefix, i.e OPDM): "  
read input_prefix  
echo "The base name is $input_prefix"  
echo  
echo "Enter inheritance pattern (AD for autosomal dominant, AR for autosomal recessive, X-linked is not available as a pre-made file): "  
read model
echo "Inheritance pattern is $model"

#model options are dominant or recessive. haven't made an x-linked file yet

perl vcf2linkdatagen.pl -annotfile annotHapMap2U.txt -pop CEU -mindepth 20 -missingness 0 -idlist $input_prefix.vcflist.txt > $input_prefix.brlmm
perl linkdatagen.pl -data m -pedfile $input_prefix.ped -whichSamplesFile $input_prefix.ws -callFile $input_prefix.brlmm -annotFile annotHapMap2U.txt -pop CEU -binsize 0.3 -mindist 0.15 -prog me -outputDir $input_prefix.HapMap2
cd $input_prefix.HapMap2_merlin/genome/
cp $model.param.tbl temp.tbl
merlin -d merlin_autosome_$input_prefix.dat -p merlin_genome_$input_prefix.ped -f merlin_genome_$input_prefix.freq -m merlin_genome_$input_prefix.map --prefix merlin_autosome_$input_prefix --error
pedwipe -d merlin_autosome_$input_prefix.dat -p merlin_genome_$input_prefix.ped -e merlin_autosome_$input_prefix.err
mv wiped.dat merlin_autosome_wiped_$input_prefix.dat
mv wiped.ped merlin_autosome_wiped_$input_prefix.ped
minx -d merlin_X_$input_prefix.dat -p merlin_genome_$input_prefix.ped -f merlin_genome_$input_prefix.freq -m merlin_genome_$input_prefix.map --prefix merlin_X_$input_prefix --error
pedwipe -d merlin_X_$input_prefix.dat -p merlin_genome_$input_prefix.ped -e merlin_X_$input_prefix.err
mv wiped.dat merlin_X_wiped_$input_prefix.dat
mv wiped.ped merlin_X_wiped_$input_prefix.ped
rm wiped.freq
mv temp.tbl $model.param.tbl
merlin -d merlin_autosome_wiped_$input_prefix.dat -p merlin_autosome_wiped_$input_prefix.ped -f merlin_genome_$input_prefix.freq -m merlin_genome_$input_prefix.map --smallswap --megabytes:9999 --founders --MarkerNames --pairs --ibd --extended --best --exp --model $model.param.tbl --pdf --tabulate --prefix merlin_autosome_$input_prefix
minx -d merlin_X_wiped_$input_prefix.dat -p merlin_X_wiped_$input_prefix.ped -f merlin_genome_$input_prefix.freq -m merlin_genome_$input_prefix.map --smallswap --megabytes:9999 --founders --MarkerNames --pairs --ibd --extended --best --exp --model $model.param.tbl --pdf --tabulate --prefix merlin_X_$input_prefix
