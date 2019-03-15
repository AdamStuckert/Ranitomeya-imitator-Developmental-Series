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

# Run the Oyster River Protocol (MacManes 2018)--this is a modified script for our readset
# Modified makefile is located in the "SupplementalDocuments" folder
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

#### Annotation with diamond to Xenopus peptide database
# download
cd ~/peptide_databases/
curl -LO ftp://ftp.ensembl.org/pub/release-95/fasta/xenopus_tropicalis/pep/Xenopus_tropicalis.JGI_4.2.pep.all.fa.gz
gunzip Xenopus_tropicalis.JGI_4.2.pep.all.fa.gz

# Make the diamond index  
diamond makedb --in Xenopus_tropicalis.JGI_4.2.pep.all.fa -d Xen_ensembl95.dmnd

# Diamond mapping
cd ~/Developmental_series/
diamond blastx -d /home/summersk/peptide_databases/Xenopus_tropicalis/Xen_ensembl95.dmnd -q Ranitomeya_imitator_transcriptome.fasta -o imi_Xen95.txt --threads 40

# sort by top hit
sort imi_Xen95.txt -k 1,1 -k11,11g | sort -u -k 1,1 --merge > imi_Xen95_tophit.txt
