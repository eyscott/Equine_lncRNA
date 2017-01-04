#**Identifying Equine lncRNA from RNA-seq data**#

The horse lacks a multi-tissue lncRNA annotation across several individuals.  We have identified candidate lncRNA from RNA-seq data
spanning eight tissues and including over fifty individuals.  We employed four levels of filtration to get these candidate lncRNA:
  1. Remove transcripts < 0.1 TPM (<5 TPM for single exon transcripts)
  2. Remove transcripts < 200 bp
  3. Remove transcripts with protein-coding capability (using hmmer and blastp)
  4. Isolate transcripts localizing to 3’& 5’ upstream sequences of genes (1 kb) 
  
All the scripts used to run this pipeline are provided in numerical order in this Github page.  The scripts used to make the figures
for the publication are also provided in this Github page. 

