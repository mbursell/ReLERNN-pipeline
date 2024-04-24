#!/bin/bash

chromlist=($(cat chroms_rj.txt))

for chrom in ${chromlist[@]}

do

count=$(bcftools view -r $chrom red_junglefowl.vcf.gz | grep -v -c '^#')

echo "$chrom:$count"

done
