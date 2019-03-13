#!/bin/bash

###################################################
###Bash script for developmental series analysis###
###################################################

# Combine all forward/reverse reads to prepare to subsample
for i in reads/*R1.fastq.gz; do zcat $i; echo; done > all_reads_R1.fastq &
for f in reads/*R2.fastq.gz; do zcat $f; echo; done > all_reads_R2.fastq 

# subsample to 40M reads for forward and reverse
seqtk sample -s100 all_reads_R1.fastq 40000000 > subsamp.R1.fastq &
seqtk sample -s100 all_reads_R2.fastq 40000000 > subsamp.R2.fastq


/home/summersk/programs/Oyster_River_Protocol/oyster35.mk main \
MEM=750 \
CPU=28 \
READ1=subsamp.R1.fastq \
READ2=subsamp.R2.fastq \
RUNOUT=subsamp

# The transcriptome is in a different file, so I'll take the 'good' transcripts and move that fasta
cp /home/summersk/Developmental_series/orthofuse/subsamp/merged/merged/good.merged.fasta .

# Rename all the transcripts, largely for aesthetics...
awk '/^>/{print ">Transcript_" ++i; next}{print}' good.merged.fasta > subsamp.imitator.merged.fasta

### Prepare individual reads
# Remove adaptors and trim the reads from each sample, first get the samples
samples=$(ls reads/*fastq.gz | sed "s/.R1.fastq//g" | sed "s/.R2.fastq//g")

# Run trimmomatic
cd reads
(ls *fastq.gz | sed "s/.R1.fastq.gz//g" | sed "s/.R2.fastq.gz//g") | \
parallel -j 10 trimmomatic-0.36.jar PE -threads 4 \
-baseout /home/summersk/Developmental_series/rcorr/{}.TRIM.fastq.gz {}.R1.fastq.gz {}.R2.fastq.gz \
LEADING:3 TRAILING:3 ILLUMINACLIP:barcodes.fa:2:30:10 MINLEN:25 

# R corrector  | sed "s/rcorr\//g"
cd ..
(ls rcorr/*P.fastq.gz | sed "s/.TRIM_1P.fastq.gz//g" | sed "s/.TRIM_2P.fastq.gz//g"  | uniq | grep -v subsamp) | \
parallel -j 6 run_rcorrector.pl -t 4 -k 31 -1 {}.TRIM_1P.fastq.gz -2 {}.TRIM_2P.fastq.gz -od rcorr


### Pseudo-quantification with Kallisto
# Make directories for each sample:
cd rcorr/
samples=$(ls *P.cor.fq.gz | sed "s/.TRIM_1P.cor.fq.gz//g" | sed "s/.TRIM_2P.cor.fq.gz//g" | uniq)
cd ..
mkdir kallisto_quants
cd kallisto_quants
for i in $samples; do mkdir $i; done
cd ..

# Pseudo-quantification with kallisto, build the index first
kallisto index -i subsamp_devseries.idx subsamp.imitator.merged.fasta

# Perform the actual pseudo-quantification

parallel -j 10 kallisto quant -i subsamp_devseries.idx -o kallisto_quants/{} -b 100 \
rcorr/{}.TRIM_1P.cor.fq.gz rcorr/{}.TRIM_2P.cor.fq.gz ::: $samples

printf " ################################################ \n \
######## Diamond annotation to Xenopus ######### \n \
################################################ \n"

#### Annotation with diamond to amphibian (xenopus, nanorana peptides) + UniRef90 databases

# cat all peptides together 
cat Nanorana_parkeri.gene.v2.pep.fa uniref90_Sept2018.fasta Xenopus_tropicalis.JGI_4.2.pep.all.fa > all_peptides.fa

# Make the diamond index  
diamond makedb --in all_peptides.fa -d allpep &

# Diamond mapping
diamond blastx -d /home/summersk/peptide_databases/allpep.dmnd -q subsamp.imitator.merged.fasta -o subimi2allpep.m8 --threads 16

# sort by top hit
sort subimi2allpep.m8 -k 1,1 -k11,11rg | sort -u -k 1,1 --merge > subimiallpep_tophit.txt

# I will also annotate just to xenopus, in case that makes a difference for pulling out gene names
diamond blastx -d /home/summersk/peptide_databases/xen.dmnd -q subsamp.imitator.merged.fasta -o subimi2xen.m8 --threads 16

# sort by top hit
sort subimi2xen.m8 -k 1,1 -k11,11g | sort -u -k 1,1 --merge > subimixen_tophit.txt


sort subimi2allpep.m8 -k 1,1 -k11,11g | sort -u -k 1,1 --merge > subimiallpep_tophit.txt


# UniRef90 hits:
grep "UniRef90_" subimiallpep_tophit.txt > tmp.txt
awk '{print $2}' tmp.txt > allpep_unirefs


# Just uniref mapping:

# Make the diamond index  
diamond makedb --in uniref90.fa -d unirefs

# Diamond mapping
diamond blastx -d /home/summersk/peptide_databases/unirefs.dmnd -q subsamp.imitator.merged.fasta -o subimi2uniref.m8 --threads 16

# Data: Total time = 2101.9s; Reported 926134 pairwise alignments, 926134 HSP; 40689 queries aligned.

# sort by top hit
sort subimi2uniref.m8 -k 1,1 -k11,11g | sort -u -k 1,1 --merge > subimiuniref_tophit.txt
# 40,689 queries aligned
# 46.1% annotation rate

# Use this to download the uniref90 hits via the website. Use IDs to get details, only taking the "reviewed" SwissProt hits!

# Now, download all of this on to my computer using cyberduck like a pleb.
	# This includes the kallisto quantifications as well as the annotation output file from diamond


## Get uniprot KB IDs from Uniref90 IDs

for i in $(cat allpep_unirefs)
do
grep $i /home/summersk/peptide_databases/uniref2uniprotkb.txt >> uniprotkbs.txt
done 


## NOTE: THIS WOULD BE WAY FASTER IN R

??Add a header to the uniref2uniprotkb file?
echo uniref_ids > tmp
cat allpep_unirefs >> tmp
mv tmp allpep_unirefs

printf "uniref_ids \t uniprotKBs \n" tmp
cat uniref2uniprotkb >> tmp


R
library(dplyr)
library(data.table)
ids <- fread("allpep_unirefsforR", header = TRUE)
kbs <- fread("/home/summersk/peptide_databases/uniref2uniprotkbforR.txt", header = TRUE)


head(ids)
head(kbs)

new <- dplyr::left_join(ids, kbs, by = "uniref_ids")

write.table(new, "uniprotkblist.tsv", sep = "\t")

q()



awk 'BEGIN{FS="\t"}{printf("%s\t%s\n", $1, $9)}' idmapping_selected.tab > NEWANNOTATION.txt
awk '{print $1}' NEWANNOTATION.txt > accessions.txt

grep "UniRef90" imi_allann_tophit.txt | awk '{printf("%s\t%s\n", $1, $2)}' > uniref_transcripts_Mar2019.txt
awk '{print $2}' uniref_transcripts_Mar2019.txt > uniref_IDS_Mar2019.txt


# query accessions.txt at https://www.uniprot.org/mapping/ with UniProtKB AC/ID to Gene name