# Data Preparation

## Create Genotype VCF

## Filter Genotype VCF
If starting whole pipeline from a publicly available gVCF that you need to filter for specific samples, begin at step 1. If your gVCF contains only the samples and chromosomes you want, start at step 3. 

Required Packages:
    * bcftools

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
5. 

## Create BED File
ReLERNN requires a BED-formatted (zero-based) file of chromosome positions for the reference genome used to create the gVCF. To create this BED file, we used the UCSC Genome Browser. These steps can be run locally on your machine. 

Required Packages:
    * bedops

1. Go to [UCSC Kent Utilities](http://hgdownload.soe.ucsc.edu/admin/exe/) and download the individual `fetchChromSizes` package into the directory of your choice using the following command
    ```
    rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/fetchChromSizes ./
    ```
2. In the same directory as `fetchChromSizes`, run this command to generate your BED file
    ```
    ./fetchChromSizes galGal4 | awk -v FS="\t" -v OFS="\t" '{ print $1, "0", $2; }' | sort-bed - > galGal4.bed
    ```
    * In this example, we generate a BED file for the galGal4 reference genome. Replace this with the name of your reference genome. This step may be limited by what reference genomes are in the UCSC Genome Browser. 

   
