# Data Preparation
This file outlines the steps we used to prepare each breed/population for ReLERNN starting from a publicly available genotype VCF. 

Packages Used:
* bcftools v1.19
* vcflib v1.0.3


1. Generated a text file of target SampleIDs
   * create a text file with one sampleID per line for each breed/population you want to extract a vcf for

2. Filtered gVCF for SampleIDs using bcftools
    ```
    bcftools view –samples-file SampleIDs.txt original.vcf > filtered.vcf
    ```
3. Zipped and indexed for next step
   ```
   bgzip filtered.vcf
   bcftools index filtered.vcf.gz
   ```
4. Filtered gVCF for chromosomes of interest using bcftools
    ```
    bcftools view filtered.vcf.gz --regions 1,2,3 > filtered_chroms.vcf
    ```
    * Note that this code selects chromosomes 1, 2, and 3 from the VCF. We changed this based on which autosomes we wanted from each species. 

5. Filtered gVCF for only biallelic variants using bcftools
    ```
    bcftools view -m2 -M2 filtered_chroms.vcf > filtered_chroms_biallelic.vcf
    ```
6. Recalculated AF and filtered out invariant sites
   ```
   vcffixup filtered_chroms_biallelic.vcf | vcffilter -f "AF > 0.05" -f "AF < 0.95" > filtered_chroms_biallelic_fixup.vcf
   ```
7. Zipped and indexed to prepare for counting snps
   ```
   bgzip filtered_chroms_biallelic_fixup.vcf
   bcftools index filtered_chroms_biallelic_fixup.vcf.gz
   ```
8. Counted the number of snps per chromosome
   * ReLERNN requires that each chromosome has at least 250 snps. Here, we count the number of snps per chromosome. If a chromosome had less that 250 snps, we removed the name of the chromosome from the BED file of chromosome positions.
   * Started by making a text file of all available chromosomes
   ```
   bcftools query -f ‘%CHROM\n’ filtered_chroms_biallelic_fixup.vcf.gz | uniq > chroms.txt
   ```
   * Ran this bash script to count snps on each chromosome in the text file.
   ```
   chmod +x snp_count.sh
   ./snp_count.sh
   ```
