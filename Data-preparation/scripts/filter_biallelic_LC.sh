#!/bin/bash

#SBATCH -p standard # partition (queue)
#SBATCH -o filter_biallelic_LC.%N.%j.out # STDOUT
#SBATCH -e filter_biallelic_LC.%N.%j.err # STDERR
#SBATCH --mem=15G

bcftools view -m2 -M2 /home6/mgbursel/BIO624/data/vcfs/local_chicken_autosomes.vcf > /home6/mgbursel/BIO624/data/vcfs/local_chicken_autosomes_biallelic.vcf
