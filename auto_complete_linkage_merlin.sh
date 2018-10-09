#! /bin/bash

file=$1
model=$2
#model options are dominant or recessive. haven't made an x-linked file yet

cd /home/san/sbeecroft/old_linkdatagen/

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
