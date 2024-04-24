#!/bin/bash

#SBATCH -p standard # partition (queue)
#SBATCH -o index_bgzip_LC.%N.%j.out # STDOUT
#SBATCH -e index_bgzip_LC.%N.%j.err # STDERR
#SBATCH --mem=10G

bcftools index /home6/mgbursel/BIO624/data/vcfs/local_chicken.vcf.gz
