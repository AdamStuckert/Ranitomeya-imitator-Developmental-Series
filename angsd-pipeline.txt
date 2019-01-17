#!/bin/bash

TRANSCRIPTOME=$"subsamp.imitator.merged.fasta"

# create an BWA index from the reference transcriptome
bwa index $TRANSCRIPTOME

# list samples to work with *these are all in the rcorr/ directory, and have already been trimmed/error corrected)
cd rcorr/
samples=$(ls *P.cor.fq.gz | sed "s/.TRIM_1P.cor.fq.gz//g" | sed "s/.TRIM_2P.cor.fq.gz//g" | uniq  | grep -v subsamp)
echo $samples

# Map individual samples to reference transcriptome using BWA
cd ..
mkdir bambams
for i in $samples; do
    bwa mem $TRANSCRIPTOME  -t 24 \
    rcorr/$i.TRIM_1P.cor.fq.gz \
    rcorr/$i.TRIM_2P.cor.fq.gz > bambams/$i.sam
done

# define picard path for next steps
PICARD=$"/home/summersk/programs/picard-2.18.5/picard.jar"

# define samples, these are the same as above
sammies=$(ls bambams/*.sam | sed "s/.sam//g")

# run picard to produce sorted sam files
parallel -j 12 java -jar $PICARD SortSam INPUT={}.sam OUTPUT={}_sorted_reads.bam SORT_ORDER=coordinate ::: $sammies

# Add group information to files
sorted=$(ls *_sorted_reads.bam | sed "s/_sorted_reads.bam//g")
parallel -j 32 java -jar $PICARD AddOrReplaceReadGroups  I={}_sorted_reads.bam O={}_marked_groups.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={} ::: $sorted

# Mark duplicates here
marked=$(ls *_marked_groups.bam | sed "s/_marked_groups.bam//g")
parallel -j 6 java -jar $PICARD MarkDuplicates INPUT={}_marked_groups.bam OUTPUT={}_dedup_reads.bam METRICS_FILE={}_metrics.txt REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 ::: $marked  


# sample list
samples=$(ls *_dedup_reads.bam)
for i in $samples
do
echo $i
done > samples4angsd.txt

### global angsd!
/home/summersk/programs/angsd/angsd -b samples4angsd.txt -anc $TRANSCRIPTOME -out angsd_global -P 20 -SNP_pval 1e-6 -minMapQ 20 -minQ 20 -setMinDepth 100 -setMaxDepth 6000 -minInd 6 -minMaf 0.01 -GL 1 -doMaf 1 -doMajorMinor 1 -doGlf 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -dosaf 1 

### Call angsd on each frog ###
# Make a SNP list with major
gunzip angsd_global.mafs.gz
cut -f 1,2,3,4 angsd_global.mafs | tail -n +2 > angsd_global_snplist.txt

# index the sites
/home/summersk/programs/angsd/angsd sites index angsd_global_snplist.txt

# run angsd on each color morph individually
cat samples4angsd.txt | grep H > Huallaga4angsd.txt
cat samples4angsd.txt | grep S > Sauce4angsd.txt
cat samples4angsd.txt | grep T > Tarapoto4angsd.txt
cat samples4angsd.txt | grep V > Varadero4angsd.txt


/home/summersk/programs/angsd/angsd -b Huallaga4angsd.txt -anc $TRANSCRIPTOME -out angsd_Huallaga_calling -P 6 -setMaxDepth 6000 -GL 1 -doMaf 1 -doMajorMinor 3 -doGlf 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -dosaf 1 -sites angsd_global_snplist.txt -minMapQ 20 -minQ 20 >& Huallaga.log &

/home/summersk/programs/angsd/angsd -b Sauce4angsd.txt -anc $TRANSCRIPTOME -out angsd_Sauce_calling -P 6 -setMaxDepth 6000 -GL 1 -doMaf 1 -doMajorMinor 3 -doGlf 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -dosaf 1 -sites angsd_global_snplist.txt -minMapQ 20 -minQ 20 >& Sauce.log &

/home/summersk/programs/angsd/angsd -b Tarapoto4angsd.txt -anc $TRANSCRIPTOME -out angsd_Tarapoto_calling -P 6 -setMaxDepth 6000 -GL 1 -doMaf 1 -doMajorMinor 3 -doGlf 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -dosaf 1 -sites angsd_global_snplist.txt -minMapQ 20 -minQ 20 >& Tarapoto.log &

/home/summersk/programs/angsd/angsd -b Varadero4angsd.txt -anc $TRANSCRIPTOME -out angsd_Varadero_calling -P 6 -setMaxDepth 6000 -GL 1 -doMaf 1 -doMajorMinor 3 -doGlf 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -dosaf 1 -sites angsd_global_snplist.txt -minMapQ 20 -minQ 20 >& Varadero.log &

