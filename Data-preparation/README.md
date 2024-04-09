# Data Preparation

## Genotype VCF Creation

## Genotype VCF Filtering
If starting whole pipeline from a publicly available gVCF that you need to filter for specific samples, begin at step 1. If your gVCF contains only the samples you want, start at step 3. 

1. Generate a text file of target SampleIDs

2. Filter gVCF for SampleIDs using bcftools
    ```
    bcftools view –samples-file SampleIDs.txt original.vcf > filtered.vcf
    ```
3. Filter gVCF for chromosomes of interest using bcftools
    ```
    bcftools view filtered.vcf --regions 1,2,3 > filtered_chroms.vcf
    ```
    * Note that this code selects chromosomes 1, 2, and 3 from the VCF. Chromsomes may be named differently in other VCFs. Make sure the names match those in the VCF. 

4. Filter gVCF for only biallelic variants using bcftools
    ```
    bcftools view -m2 -M2 filtered_chroms.vcf > filtered_chroms_biallelic.vcf
    ```
