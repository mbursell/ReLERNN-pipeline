#!/bin/bash

#SBATCH -p standard # partition (queue)
#SBATCH -o extract_YC.%N.%j.out # STDOUT
#SBATCH -e extract_YC.%N.%j.err # STDERR
#SBATCH --exclusive
#SBATCH --mem=0

bcftools view --samples-file /home6/mgbursel/BIO624/data/sampleIDs/yunbao_chicken_samples.txt /home6/mgbursel/download.big.ac.cn/chickensd/vcf/gz/Global_chicken.GA_without.merge.raw.filter4_final_no3.vcf > /home6/mgbursel/BIO624/data/vcfs/yunbao_chicken.vcf 
