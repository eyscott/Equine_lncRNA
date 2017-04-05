**Identifying Equine lncRNA from RNA-seq data**

The horse lacks a multi-tissue lncRNA annotation across several individuals.  We have identified candidate lncRNA from RNA-seq data
spanning eight tissues and including over fifty individuals.  We employed four levels of filtration to get these candidate lncRNA:
  1. Remove transcripts < 0.1 TPM (<5 TPM for single exon transcripts)
  2. Remove transcripts < 200 bp
  3. Remove transcripts with protein-coding capability (using hmmer and blastp)
  4. Isolate transcripts localizing to 3’& 5’ upstream sequences of genes (1 kb) 
  
All the scripts used to run this pipeline are provided in numerical order in this Github page.  The scripts used to make the figures
for the publication are also provided, ordered using roman numerals, in this Github page. 

The input data originated from the equine transcriptome annotation by Mansour et al (2016). Scripts and data needed for the inputs can be found 
[here](https://github.com/drtamermansour/horse_trans) 

The lncRNA gene transfer format (GTF) file can be downloaded [here] (http://de.cyverse.org/dl/d/EBBAF4AF-7A5D-42F2-AE59-E1FF5A3C460C/lncRNA_new.gtf)

The lncRNA browser extensible data (BED) file can be downloaded [here] (http://de.cyverse.org/dl/d/7A350F1A-C796-4D0C-B2DF-36AC604BBE11/lncRNA_final.bed)  

For a full equine transcriptome GTF file where the lncRNA and protein-coding transcripts are labeled in the "gene_type" section (which is column 9 of the GTF) either as "non_coding" or "protein_coding" [this] (http://de.cyverse.org/dl/d/5EA00FE4-4581-4A06-A6AB-1103AD56F6F0/ecTranscriptome.gtf) file

The BED file can be added as a custom track in UCSC browser for visualization

