cd ~/Desktop/lncRNA
wget ftp://ftp.ensembl.org/pub/release-80/gtf/equus_caballus/Equus_caballus.EquCab2.80.gtf.gz
gunzip Equus_caballus.EquCab2.80.gtf.gz
grep -v "^#" Equus_caballus.EquCab2.80.gtf | awk -F "\t" '{if($3=="transcript")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "\""); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 2)] = substr($i, n + 1, length(x) - 1) } } gene_id = vars["gene_id"]; gene_version = vars[" gene_version"]; transcript_id = vars[" transcript_id"]; transcript_version = vars[" transcript_version"]; gene_biotype = vars[" gene_biotype"]; } { print gene_id "." gene_version, transcript_id "." transcript_version, gene_biotype }' | awk -F"\t" '{if($3=="misc_RNA")print}'> ENS.ncRNA;
awk -F"\t" '{print $2}' ENS.ncRNA | grep -Fwf - inputs/annotations/RNAseqSupTrans.merge.reduced > RNAseqSupTrans.merge.reduced.ENS.ncRNA

wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/GFF/ref_EquCab2.0_top_level.gff3.gz
gunzip ref_EquCab2.0_top_level.gff3.gz
grep -v "^#" ref_EquCab2.0_top_level.gff3 | awk -F "\t" '{if($3=="ncRNA")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } gene = vars["gene"]; ncrna_class = vars["ncrna_class"]; transcript_id = vars["transcript_id"]; } { print gene,transcript_id,ncrna_class }' | awk -F"\t" '{if($3=="lncRNA")print}' > NCBI.ncRNA;
awk -F"\t" '{print $2}' NCBI.ncRNA | grep -Fwf - inputs/annotations/RNAseqSupTrans.merge.reduced > RNAseqSupTrans.merge.reduced.NCBI.ncRNA

head -n1 inputs/annotations/RNAseqSupTrans.merge.reduced > RNAseqSupTrans.merge.reduced.ncRNA
cat RNAseqSupTrans.merge.reduced.ENS.ncRNA RNAseqSupTrans.merge.reduced.NCBI.ncRNA | sort | uniq >> RNAseqSupTrans.merge.reduced.ncRNA


