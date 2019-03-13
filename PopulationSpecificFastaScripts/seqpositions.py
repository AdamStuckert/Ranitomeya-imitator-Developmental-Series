#!/usr/bin/env python3

######
# purpose: 
# input: fasta file, list of positions to change, list of nucleotides to change to, name of output file
# output: fasta file
# return: fasta file with different nucleotides at positions
# name of this python script: seqpositions.py
######

import csv

with open('imitator_fixed_SNPs_colorgenes.csv', mode='r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    line_count = 0
    for row in csv_reader:
        gene = {row["full_gene_name"]}
        contig = {row["contig"]}
        position = {row["position"]}
        if {row["Huallaga_minor_frequency"]} == 1:
            Huallaga_nuc = {row["Huallaga_minor"]}
        else:
            Huallaga_nuc = {row["Huallaga_major"]}
        if {row["Sauce_minor_frequency"]} == 1:
            Sauce_nuc = {row["Sauce_minor"]}
        else:
            Sauce_nuc = {row["Sauce_major"]}
        if {row["Tarapoto_minor_frequency"]} == 1:
            Tarapoto_nuc = {row["Tarapoto_minor"]}
        else:
            Tarapoto_nuc = {row["Tarapoto_major"]}
        if {row["Varadero_minor_frequency"]} == 1:
            Varadero_nuc = {row["Varadero_minor"]}
        else:
            Varadero_nuc = {row["Varadero_major"]}

        # clean up the language    
        gene = gene.replace("{'", "")
        gene = gene.replace("'}", "") 
        position = position.replace("{'", "")
        position = position.replace("'}", "")
        Huallaga_nuc = Huallaga_nuc.replace("{'", "")
        Huallaga_nuc = Huallaga_nuc.replace("'}", "")
        Sauce_nuc = Sauce_nuc.replace("{'", "")
        Sauce_nuc = Sauce_nuc.replace("'}", "")
        Tarapoto_nuc = Tarapoto_nuc.replace("{'", "")
        Tarapoto_nuc = Tarapoto_nuc.replace("'}", "")
        Varadero_nuc = Varadero_nuc.replace("{'", "")
        Varadero_nuc = Varadero_nuc.replace("'}", "")
        print("{0} \t{1} \t{2} \t{3} \t{4} \t{5} \t{6}".format(gene, position, Huallaga_nuc, Sauce_nuc, Tarapoto_nuc, Varadero_nuc))
        line_count += 1

    print(f'Processed {line_count} lines.')