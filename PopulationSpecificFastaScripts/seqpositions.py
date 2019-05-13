#!/usr/bin/env python3

######
# purpose: 
# input: fasta file, list of positions to change, list of nucleotides to change to, name of output file
# output: fasta file
# return: fasta file with different nucleotides at positions
######

import csv

with open('imitator_fixed_SNPs_colorgenes.csv', mode='r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    line_count = 0
    print("Gene \tPosition \tHuallaga \tSauce \tTarapoto \tVaradero")
    for row in csv_reader:
        gene = str(row["gene_name"])
        contig = str(row["contig"])
        position = str(row["position"])
        if row["Huallaga_minor_frequency"] == "1":
            Huallaga_nuc = str(row["Huallaga_minor"])
        else:
            Huallaga_nuc = str(row["Huallaga_major"])
        if row["Sauce_minor_frequency"] == "1":
            Sauce_nuc = str(row["Sauce_minor"])
        else:
            Sauce_nuc = str(row["Sauce_major"])
        if row["Tarapoto_minor_frequency"] == "1":
            Tarapoto_nuc = str(row["Tarapoto_minor"])
        else:
            Tarapoto_nuc = str(row["Tarapoto_major"])
        if row["Varadero_minor_frequency"] == "1":
            Varadero_nuc = str(row["Varadero_minor"])
        else:
            Varadero_nuc = str(row["Varadero_major"])

        print("{0} \t{1} \t{2} \t{3} \t{4} \t{5} \t{6}".format(gene, position, Huallaga_nuc, Sauce_nuc, Tarapoto_nuc, Varadero_nuc, contig))
        line_count += 1

