#!/bin/bash

echo 'Starting downloads. This will take a few minutes'
wget http://bioinf.wehi.edu.au/software/linkdatagen/linkdatagen.pl
wget http://bioinf.wehi.edu.au/software/linkdatagen/vcf2linkdatagen.pl
wget http://bioinf.wehi.edu.au/software/linkdatagen/annotation/affymetrix/mappingfiles_affy_build37.tgz
wget http://bioinf.wehi.edu.au/software/linkdatagen/annotation/annotIllumina.tgz
wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
wget http://bioinf.wehi.edu.au/software/linkdatagen/annotation/mps/annotHapMap3U.txt.gz
wget http://bioinf.wehi.edu.au/software/linkdatagen/annotation/mps/annotHapMap2U.txt.gz
conda install -c bioconda merlin 

echo 'Decompressing files'
tar -xzvf mappingfiles_affy_build37.tgz
tar -xvzf annotIllumina.tgz
tar -xzvf merlin-1.1.2.tar.gz
tar -xf samtools-0.1.19.tar.bz2
gunzip annotHapMap2U.txt.gz 
gunzip annotHapMap3U.txt.gz
gunzip hg19.fa.gz

echo '''
==============READ ME===================
Should be complete, check the outputs

To install Merlin using simple defaults, do the following:

cd merlin-1.1.2
make 
make install

=======================================

To install samtools using simple defaults, do the following:

cd samtools-0.1.19
make 
'''

