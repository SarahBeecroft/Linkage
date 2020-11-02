#!/bin/bash

for number in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
do 

cut -f2,6-12 maori_homozygosity_mapping.HapMap2_complete/cp_$number.ped >> maori_homozygosity_mapping.HapMap2_complete/cp_ALL_edited.ped

done

