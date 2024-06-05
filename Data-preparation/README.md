# Data Preparation
ReLERNN requires a genotype VCF that only contains biallelic variants and a BED file of the reference genome chromosome positions. Here we will show the steps we used starting from a publicly available gVCF. 

## Subset Samples from Genotype VCF
If you are beginning with a Genotype VCF that contains more than your target samples, begin at this step. If you want to include all samples in your VCF, proceed to the next section. 

Required Packages:
   * bcftools (we installed via conda and mamba, but feel free to use another method if you prefer another)

1. Generate a text file of target SampleIDs
   * create a text file with one sampleID per line for each breed/population you want to extract a vcf for
   
2. Filter gVCF for SampleIDs using bcftools
    ```
    bcftools view –samples-file SampleIDs.txt original.vcf > filtered.vcf
    ```

## Filter Genotype VCF
If starting from an available gVCF that you need to filter for specific chromosomes, begin at step 1. If your gVCF already contains only the chromosomes you want, start at step 3. 

Required Packages:
   * bcftools (we installed via conda and mamba, but feel free to use another method if you prefer another)
   * vcflib v1.0.3

1. Zip and index for next step
   ```
   bgzip filtered.vcf
   bcftools index filtered.vcf.gz
   ```
2. Filter gVCF for chromosomes of interest using bcftools
    ```
    bcftools view filtered.vcf.gz --regions 1,2,3 > filtered_chroms.vcf
    ```
    * Note that this code selects chromosomes 1, 2, and 3 from the VCF. Chromsomes may be named differently in other VCFs. Make sure the names match those in the VCF. 

3. Filter gVCF for only biallelic variants using bcftools
    ```
    bcftools view -m2 -M2 filtered_chroms.vcf > filtered_chroms_biallelic.vcf
    ```
4. Recalculate AF and filter out invariant sites
   ```
   vcffixup filtered_chroms_biallelic.vcf | vcffilter -f "AF > 0.05" -f "AF < 0.95" > filtered_chroms_biallelic_fixup.vcf
   ```
5. Zip and index to prepare for counting snps
   ```
   bgzip filtered_chroms_biallelic_fixup.vcf
   bcftools index filtered_chroms_biallelic_fixup.vcf.gz
   ```
6. Count the number of snps per chromosome
   * ReLERNN requires that each chromosome has at least 250 snps. Here, we count the number of snps per chromosome. If a chromosome has less that 250 snps, you will need to remove the name of the chromosome from the BED file we generate in the section below.
   * Start by making a text file of all available chromosomes
   ```
   bcftools query -f ‘%CHROM\n’ filtered_chroms_biallelic_fixup.vcf.gz | uniq > chroms.txt
   ```
   * Run this bash script to count snps on each chromosome in the text file. Make sure to edit the name of the text file and vcf in the bash script
   ```
   chmod +x snp_count.sh
   ./snp_count.sh
   ```
    
## Create BED File
ReLERNN requires a BED-formatted (zero-based) file of chromosome positions for the reference genome used to create the gVCF. To create this BED file, we used the UCSC Genome Browser. These steps can be run locally on your machine. 

Required Packages (tested with conda and mamba installation method):
   * bedops
   * python

1. Go to [UCSC Kent Utilities](http://hgdownload.soe.ucsc.edu/admin/exe/) and download the individual `fetchChromSizes` package into the directory of your choice using the following command
    ```
    rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/fetchChromSizes ./
    ```
2. In the same directory as `fetchChromSizes`, run this command to generate your BED file
    ```
    ./fetchChromSizes galGal4 | awk -v FS="\t" -v OFS="\t" '{ print $1, "0", $2; }' | sort-bed - > galGal4.bed
    ```
    * In this example, we generate a BED file for the galGal4 reference genome. Replace this with the name of your reference genome. This step may be limited by what reference genomes are in the UCSC Genome Browser.
   
3. Filter BED file for only chromosomes of interest
   * The initial BED file will contain all chromosomes in the reference, including sex chromosomes and microchromosomes. You will need to filter for only your chromosomes of interest.
   * This step is highly dependent on what data you are working with. For the chickens, the BED file used the formatting "chr1" and the VCF used "1" to label for chromosome 1. Smaller chromosomes had additional labeling after "chr" and the chromosome number. We used a python script to only grab the large chromosomes with no additional labelling using regular expressions.
   * Warning: The names of the chromosomes in the BED file and the VCF could be different. Make sure they are matching.
   ```
   python Filter_BED.py
   ```
   * If any of the chromosomes had less than 250 snps, please remove the name of that chromosome from the BED file before proceeding to run ReLERNN.
  
4. Filter BED file for structural variants

   
