#!/bin/bash


## This is called poplists.sh

./seqpositions.py  > position_list.txt

awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' position_list.txt > Huallaga_SNP_list.txt
awk 'BEGIN{OFS="\t"} {print $1,$2,$4}' position_list.txt > Sauce_SNP_list.txt
awk 'BEGIN{OFS="\t"} {print $1,$2,$5}' position_list.txt > Tarapoto_SNP_list.txt
awk 'BEGIN{OFS="\t"} {print $1,$2,$6}' position_list.txt > Varadero_SNP_list.txt

pops=$(ls *SNP_list.txt)
genes=$(awk '{print $1}' position_list.txt | sort | uniq)

# print positions for analysis by gene
for gene in $genes
do
grep $gene position_list.txt | awk '{print $2}' > ${gene}_positions.txt
done

# write nucleotide positions to a new file
for gene in $genes
do
for pop in $pops
do
grep $gene $pop | awk '{print $3}' > ${gene}_${pop}_nucs.txt
done
done

# write contig for each gene to a new fasta file
for gene in $genes
do
contig=$(grep $gene position_list.txt | awk '{print $7}' | sort | uniq)
nucs=$(sed -n "/${contig}/,/>Transcript/{/>Transcript/d;p}" /home/summersk/Developmental_series/Ranitomeya_imitator_transcriptome.fasta) 
printf '>%s \n%s\n' "$gene" "$nucs" > ${gene}.fa
done

# do all the actual reassemblies
# write a quick wrapper function
function DoTheNucSub() {
gene=$1
./rewrite_seqs.py ${gene}.fa Huallaga ${gene}_positions.txt ${gene}_Huallaga_SNP_list.txt_nucs.txt ${gene}_Huallaga.fa

./rewrite_seqs.py ${gene}.fa Sauce ${gene}_positions.txt ${gene}_Sauce_SNP_list.txt_nucs.txt ${gene}_Sauce.fa

./rewrite_seqs.py ${gene}.fa Tarapoto ${gene}_positions.txt ${gene}_Tarapoto_SNP_list.txt_nucs.txt ${gene}_Tarapoto.fa

./rewrite_seqs.py ${gene}.fa Varadero ${gene}_positions.txt ${gene}_Varadero_SNP_list.txt_nucs.txt ${gene}_Varadero.fa

cat ${gene}.fa ${gene}_Huallaga.fa ${gene}_Sauce.fa ${gene}_Tarapoto.fa ${gene}_Varadero.fa > combined_${gene}.fa
}


for gene in $gene
do
DoTheNucSub $gene
done