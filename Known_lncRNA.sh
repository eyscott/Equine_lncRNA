cd ~/Desktop/lncRNA
wget ftp://ftp.ensembl.org/pub/release-80/gtf/equus_caballus/Equus_caballus.EquCab2.80.gtf.gz
gunzip Equus_caballus.EquCab2.80.gtf.gz
grep -v "^#" Equus_caballus.EquCab2.80.gtf | awk -F "\t" '{if($3=="transcript")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "\""); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 2)] = substr($i, n + 1, length(x) - 1) } } gene_id = vars["gene_id"]; gene_version = vars[" gene_version"]; transcript_id = vars[" transcript_id"]; transcript_version = vars[" transcript_version"]; gene_biotype = vars[" gene_biotype"]; } { print gene_id "." gene_version, transcript_id "." transcript_version, gene_biotype }' | awk -F"\t" '{if($3=="misc_RNA")print}'> ENS.ncRNA;
awk -F"\t" '{print $2}' ENS.ncRNA | sort | uniq | wc -l ## 175
awk -F"\t" '{print $2}' ENS.ncRNA | grep -Fwf - inputs/annotations/RNAseqSupTrans.merge.reduced > RNAseqSupTrans.merge.reduced.ENS.ncRNA ## 33
cat RNAseqSupTrans.merge.reduced.ENS.ncRNA | awk -F"\t" '{print $19}' | sort | uniq | wc -l  ## 18

wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/GFF/ref_EquCab2.0_top_level.gff3.gz
gunzip ref_EquCab2.0_top_level.gff3.gz
grep -v "^#" ref_EquCab2.0_top_level.gff3 | awk -F "\t" '{if($3=="ncRNA")print $9}' | awk -F ";" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } gene = vars["gene"]; ncrna_class = vars["ncrna_class"]; transcript_id = vars["transcript_id"]; } { print gene,transcript_id,ncrna_class }' | awk -F"\t" '{if($3=="lncRNA")print}' > NCBI.ncRNA;
awk -F"\t" '{print $2}' NCBI.ncRNA | sort | uniq | wc -l ## 4369
awk -F"\t" '{print $2}' NCBI.ncRNA | grep -Fwf - inputs/annotations/RNAseqSupTrans.merge.reduced > RNAseqSupTrans.merge.reduced.NCBI.ncRNA ## 3932
cat RNAseqSupTrans.merge.reduced.NCBI.ncRNA | awk -F"\t" '{print $15}' | sort | uniq | wc -l  ## 1720

head -n1 inputs/annotations/RNAseqSupTrans.merge.reduced > RNAseqSupTrans.merge.reduced.ncRNA
cat RNAseqSupTrans.merge.reduced.ENS.ncRNA RNAseqSupTrans.merge.reduced.NCBI.ncRNA | sort | uniq >> RNAseqSupTrans.merge.reduced.ncRNA ## 3956

tail -n+2 RNAseqSupTrans.merge.reduced.ncRNA | awk -F"\t" '{print $19}' | sort | uniq > ENS.ids ## 476
tail -n+2 RNAseqSupTrans.merge.reduced.ncRNA | awk -F"\t" '{print $15}' | sort | uniq > NCBI.ids ## 1726
comm -12 ENS.ids <(awk -F"\t" '{print $2}' ENS.ncRNA | sort | uniq) | wc -l ## 18 (the rest of the 476 ids are coding ENS transcripts assigned)
comm -12 NCBI.ids <(awk -F"\t" '{print $2}' NCBI.ncRNA | sort | uniq) | wc -l ## 1720 (the rest of 1726 ids are coding NCBI transcripts)

tail -n+2 RNAseqSupTrans.merge.reduced.ncRNA | awk -F"\t" '{print $1}' | sort > refined_ncRNA.ids ## 3956
cat inputs/allTissues_BED/supported.bed | awk -F"\t" '{print $4}' | sort > ALLsupported.ids
cat inputs/allTissues_BED/refined.bed | awk -F"\t" '{print $4}' | sort > ALLrefined.ids
comm -12 refined_ncRNA.ids ALLsupported.ids | wc -l ## 3956
comm -12 refined_ncRNA.ids ALLrefined.ids | wc -l  ## 3956

tail -n+2 inputs/novelAnn/sup/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.sup | awk -F"\t" '{print $1}' | sort > novel_I.ids 
tail -n+2 inputs/novelAnn/unsup.cons/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.cons | awk -F"\t" '{print $1}' | sort > novel_II.ids
tail -n+2 inputs/novelAnn/unsup.uncons.ORF/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.uncons.ORF | awk -F"\t" '{print $1}' | sort > novel_III.ids
comm -12 refined_ncRNA.ids novel_I.ids > novel_I_ncRNA.ids ## 2634
comm -12 refined_ncRNA.ids novel_II.ids > novel_II_ncRNA.ids ## 117
comm -12 refined_ncRNA.ids novel_III.ids > novel_III_ncRNA.ids ## 136
comm -12 refined_ncRNA.ids inputs/intergenic/unsup.ids > intergenic_ncRNA.ids ## 0
######
Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T,sep="\t");'\
'candidate_novels=subset(data1,NCBI.class_code %in% "u" | ensGTF_file.class_code %in% "u");'\
'write.table(candidate_novels,"RNAseqSupTrans.merge.reduced.candidate_novels.ncRNA", sep="\t", quote=F, row.names=F, col.names=T);' RNAseqSupTrans.merge.reduced.ncRNA
tail -n+2 RNAseqSupTrans.merge.reduced.candidate_novels.ncRNA | awk -F"\t" '{print $1}' | sort > candidate_novels_ncRNA.ids ## 2887
comm -23 refined_ncRNA.ids candidate_novels_ncRNA.ids > known_ncRNA.ids

#Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T,sep="\t");'\
#'check=c("=","j","o","x","c");'\
#'supporting_novels=subset(data1,NCBI.class_code %in% check | ensGTF_file.class_code %in% check | ISME.PBMC.class_code %in% check | Hestand_2014.class_code %in% check);'\
#'write.table(supporting_novels,"RNAseqSupTrans.merge.reduced.novel_I.ncRNA", sep="\t", quote=F, row.names=F, col.names=T);' RNAseqSupTrans.merge.reduced.candidate_novels.ncRNA
#tail -n+2 RNAseqSupTrans.merge.reduced.novel_I.ncRNA | awk -F"\t" '{print $1}' | sort > novel_I_ncRNA.ids ## 2634


#Rscript -e 'args=(commandArgs(TRUE)); data1=read.table(args[1],header=T,sep="\t");'\
#'check=c("=","j","o","x","c");'\
#'bothAnn=subset(data1,(NCBI.class_code %in% check & ensGTF_file.class_code %in% check));'\
#'write.table(bothAnn,"RNAseqSupTrans.merge.reduced.both.ncRNA", sep="\t", quote=F, row.names=F, col.names=T);' RNAseqSupTrans.merge.reduced.ncRNA
#tail -n+2 RNAseqSupTrans.merge.reduced.both.ncRNA | awk -F"\t" '{print $1}' | sort > bothAnn_ncRNA.ids ## 920

