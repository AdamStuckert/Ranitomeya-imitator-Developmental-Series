# *Ranitomeya imitator* developmental series 

In this project, we attempt to examine genes that contribute to the production of color and pattern phenotypes in the mimic poison frog *Ranitomeya imitator.* We use RNA seq data collected from tadpoles at various stages of development from four different color morphs. 

Transcriptome assembly, read cleaning, and transcript expression quantification is all done in bash. Subsequent differential gene expression analyses are in the R markdown files.

## General workflow

We assembled a *de novo* transcriptome and annotated it. Subsequently we cleaned and trimmed individual reads and mapped each sample to the assembled transcriptome. After this, we conducted differential expression analyses to examine changes in gene expression over time and between color morphs. Additionally, we examined genomic variants by calling SNPs, and looked at fixed differences between color morphs. Finally, we conducted statistical overrepresentation analyses in Panther (no code for that here, but we have included the list of genes we loaded into Panther). Finally, we produced figures and tables.

### Script order

1. 


#### Supplementary documents 

Documents not produced by these scripts are included in a separate supplemental documents folder. This includes our assembled transcriptome(?), sample metadata, etc.
