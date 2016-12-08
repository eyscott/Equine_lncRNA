module load GNU/4.4.5
module load BEDTools/2.24.0
##get fasta files for our lncRNA
bedtools getfasta -s -fi equCab2.0_genome.fa -bed lncRNA_final.bed -name -split -fo lncRNA_final.fa
##then our transcriptome
bedtools getfasta -s -fi equCab2.0_genome.fa -bed refined_nolncRNA.bed -name -split -fo refined_nolncRNA.fa
#then all our associated promoters
awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3}' lncRNA_promoters.bed | sort | uniq > lncRNA_promoters_uniq.bed
bedtools getfasta -s -fi equCab2.0_genome.fa -bed lncRNA_promoters_uniq.bed -fo lncRNA_promoters.fa

awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3}' lncRNA_down.bed | sort | uniq > lncRNA_down_uniq.bed
bedtools getfasta -s -fi equCab2.0_genome.fa -bed lncRNA_down_uniq.bed -fo lncRNA_down.fa

awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3}' refined_nolncRNA_5.bed | sort | uniq > refined_nolncRNA_5_uniq.bed
bedtools getfasta -s -fi equCab2.0_genome.fa -bed refined_nolncRNA_5_uniq.bed -fo refined_nolncRNA_5.fa

awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3}' refined_nolncRNA_3.bed | sort | uniq > refined_nolncRNA_3_uniq.bed
bedtools getfasta -s -fi equCab2.0_genome.fa -bed refined_nolncRNA_3_uniq.bed -fo refined_nolncRNA_3.fa

#Make blast db for human ncRNA,genes and genome
module load GNU/4.8.3
module load BLAST+/2.3.0
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gunzip Homo_sapiens.GRCh38.cds.all.fa.gz
cat Homo_sapiens.GRCh38.ncrna.fa Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens.GRCh38.trans.fa
makeblastdb -in Homo_sapiens.GRCh38.trans.fa -input_type fasta -dbtype nucl -out hg38_trans_db
blastn -query lncRNA_final.fa  -db hg38_trans_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > lncRNA_hg_trans.outfmt6
blastn -query refined_nolncRNA.fa  -db hg38_trans_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > coding_hg_trans.outfmt6

wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz
makeblastdb -in Homo_sapiens.GRCh38.dna_rm.toplevel.fa -input_type fasta -dbtype nucl -out hg38_genome_db


#makeblastdb -in Homo_sapiens.GRCh38.ncrna.fa -input_type fasta -dbtype nucl -out hg38_db
#makeblastdb -in Homo_sapiens.GRCh38.cds.all.fa -input_type fasta -dbtype nucl -out hg38_genes_db
#makeblastdb -in Homo_sapiens.GRCh38.dna.alt.fa -input_type fasta -dbtype nucl -out hg38_genome_db

#blast our genes and lncRNA to respective human databases
#blastn -query lncRNA_final.fa  -db hg38_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > lncRNA_hg_nt.outfmt6
#blastn -query refined_nolncRNA.fa  -db hg38_genes_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > human_genes_nt.outfmt6

#blastn -query lncRNA_final.fa  -db hg38_genome_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > lncRNA_humanGenome_nt.outfmt6
#blastn -query refined_nolncRNA.fa  -db hg38_genome_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > genes_humanGenome_nt.outfmt6

#do the same for our promoters
blastn -query lncRNA_promoters.fa  -db hg38_genome_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > lncRNApromoter_human_nt.outfmt6
blastn -query lncRNA_down.fa  -db hg38_genome_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > lncRNAdown_human_nt.outfmt6
blastn -query refined_nolncRNA_5.fa  -db hg38_genome_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > genepromoter_human_nt.outfmt6
blastn -query refined_nolncRNA_3.fa  -db hg38_genome_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > genedown_human_nt.outfmt6

