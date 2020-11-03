#!/bin/bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/twoBitToFa
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit
chmod 777 twoBitToFa
./twoBitToFa hg19.2bit ucsc.hg19.fasta
samtools faidx ucsc.hg19.fasta
bwa index -b 10gb ucsc.hg19.fasta

gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.amb
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.ann
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.bwt
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.pac
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.sa
