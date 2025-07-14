#vcf read and parse

import os
import gzip
import re
import numpy as np

def read_vcf_file(vcf_file_path):
    
    #determine the file extension
    file_extension = os.path.splitext(vcf_file_path)[1]
    
    if file_extension == ".gz":
        with gzip.open(vcf_file_path, "rt") as vcf:
            return process_vcf_file(vcf)
    elif file_extension == ".vcf":
        with open(vcf_file_path, "r") as vcf:
            return process_vcf_file(vcf)
    else:
        print("Unsupported file type")
        return [], []


def process_vcf_file(file):
    #header_printed = False
    #head = []
    #map_list = []
    position_list = []
    genotype_list = []

    for line in file:
        # Skip header lines (lines starting with '##')
        if line.startswith("#"):
            continue
            #head.append(line)
        # Process the header line (line starting with '#')
        #elif line.startswith("#"):
        #   headers = line.strip().split("\t")
        #   sample_names = headers[9:]  # Sample columns start from index 9
        #   header = ["taxa"] + sample_names
        #   #vcf_data.append("\t".join(header))
        #   header_printed = True
        #   continue
        else:
            # Split the line into columns
            columns = line.strip().split("\t")
            chrom = columns[0]
            pos = columns[1]
            genotypes = columns[9:]  # Genotype information starts from index 9
            position_list.append([chrom, pos])
            genotype_list.append(genotypes)
    return position_list, genotype_list

# Define a function to process each genotype according to the given rules
def replacement (genotype):
    
    if '.' in genotype:
        #return 'NA'
        return np.nan
    
    # use regex to split by either "|" or "/"
    parts = re.split(r'[|/]', genotype)
    #ensure there are exactly two parts after splitting.
    if len(parts) == 2:
        left, right = parts
        if left == right:
            if left == '0':
                return "0"
            else:
                return "2"
        else:
            return "1"
        
    return genotype
    

# Function to process the input file
def process_genotype(genotype_data):
    #with open(genotype_data, "r") as genotype:
    processed_data = []
    #for line in genotype_data:
    for line in genotype_data:
        genotypes = [replacement(g) for g in line]
        processed_data.append(genotypes)
        
    return processed_data
