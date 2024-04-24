#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:36:44 2024
This script takes the output BED file from UCSC Genome Browser and cleans/filters it. We remove random
chromosomes, sex chromosomes, and shorter microchromosomes. We also rename the chromosomes so they match the VCF file. 

@author: maddybursell
"""
            
            
import re

# Define the output file name
output_file = "galGal4_filtered.bed"

# Open the BED file for reading and the output file for writing
with open("galGal4.bed", "r") as bed_file, open(output_file, "w") as output:
    # Iterate through each line in the file
    for line in bed_file:
        # Split the line by tabs to get individual columns
        columns = line.strip().split("\t")
        # Extract the name from the first column
        name = columns[0]
        # Check if the name matches "chr" followed by a number using regular expression
        if re.match(r'^chr(\d+)$', name):
            # Extract the number part
            number = re.match(r'^chr(\d+)$', name).group(1)
            # Rename the name by removing "chr"
            new_name = number
            # Modify the first column with the new name
            columns[0] = new_name
            # Join the columns back together with tabs
            modified_line = "\t".join(columns)
            # Write the modified line to the output file
            output.write(modified_line + "\n")


