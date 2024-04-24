#!/bin/bash

#SBATCH -p standard # partition (queue)
#SBATCH -o extract_autosomes_YC.%N.%j.out # STDOUT
#SBATCH -e extract_autosomes_YC.%N.%j.err # STDERR
#SBATCH --mem=20G

bcftools view /home6/mgbursel/BIO624/data/vcfs/yunbao_chicken.vcf.gz --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,32 > /home6/mgbursel/BIO624/data/vcfs/yunbao_chicken_autosomes.vcf
