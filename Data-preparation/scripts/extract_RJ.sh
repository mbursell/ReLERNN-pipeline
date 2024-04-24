#!/bin/bash

#SBATCH -p standard # partition (queue)
#SBATCH -o extract_RJ.%N.%j.out # STDOUT
#SBATCH -e extract_RJ.%N.%j.err # STDERR
#SBATCH --exclusive
#SBATCH --mem=0

bcftools view --samples-file /home6/mgbursel/BIO624/data/sampleIDs/red_junglefowl_samples.txt /home6/mgbursel/download.big.ac.cn/chickensd/vcf/gz/Global_chicken.GA_without.merge.raw.filter4_final_no3.vcf > /home6/mgbursel/BIO624/data/vcfs/red_junglefowl.vcf 
