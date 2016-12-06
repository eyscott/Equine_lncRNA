#step 1: download Pfam-A db from:  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

#step 2: get the fasta format of the lncRNA from the .bed file:
#using bedtools to get the fasta
##main script for all novel categories passing the filters so far
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab2/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
tar xvzf chromFa.tar.gz
cat *.fa > equCab2.0_genome.fa
module load BEDTools/2.24.0  
bedtools getfasta -s -fi equCab2.0_genome.fa -bed novel_I_f2.bed -fo novel_I_f2.fa
bedtools getfasta -s -fi equCab2.0_genome.fa -bed novel_II_f2.bed -fo novel_II_f2.fa
bedtools getfasta -s -fi equCab2.0_genome.fa -bed novel_III_f2.bed -fo novel_III_f2.fa
bedtools getfasta -s -fi equCab2.0_genome.fa -bed intergenic_f2.bed -fo intergenic_f2.fa
bedtools getfasta -s -fi equCab2.0_genome.fa -bed known_ncRNA_f2.bed -fo known_ncRNA_f2.fa

#step 3: get the protein coding/ORF sequences of your data using Transdecoder (this takes quite some time):
#main
module load TransDecoder/3.0.0
TransDecoder.LongOrfs -t novel_I_f2.fa -m 20
TransDecoder.LongOrfs -t novel_II_f2.fa -m 20
TransDecoder.LongOrfs -t novel_III_f2.fa -m 20
TransDecoder.LongOrfs -t intergenic_f2.fa -m 20
TransDecoder.LongOrfs -t known_ncRNA_f2.fa -m 20 
#with -m 20, meaning ORF minimum of 20 (default=100) because we are dealing with lncRNA

##step 4: Predict likely coding regions with HMMER: hmmsearch, using ORFs predicted above:
module load HMMER/3.1b2
hmmsearch --cpu 8 --tblout novel_I_pfam_new.tblout Pfam-A.hmm novel_I_f2.fa.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu 8 --tblout novel_II_pfam_new.tblout Pfam-A.hmm novel_II_f2.fa.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu 8 --tblout novel_III_pfam_new.tblout Pfam-A.hmm novel_III_f2.fa.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu 8 --tblout intergenic_pfam_new.tblout Pfam-A.hmm intergenic_f2.fa.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu 8 --tblout known_lncRNA_pfam_new.tblout Pfam-A.hmm known_ncRNA_f2.fa.transdecoder_dir/longest_orfs.pep

###step 5 running blastp
module load GNU/4.8.3
module load BLAST+/2.3.0
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot
blastp -query novel_I_f2.fa.transdecoder_dir/longest_orfs.pep  -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > novel_I_sprot.outfmt6
blastp -query novel_II_f2.fa.transdecoder_dir/longest_orfs.pep  -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > novel_II_sprot.outfmt6
blastp -query novel_III_f2.fa.transdecoder_dir/longest_orfs.pep  -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > novel_III_sprot.outfmt6
blastp -query intergenic_f2.fa.transdecoder_dir/longest_orfs.pep  -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > intergenic_sprot.outfmt6
blastp -query known_ncRNA_f2.fa.transdecoder_dir/longest_orfs.pep  -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > known_lncRNA_sprot.outfmt6

wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

grep "gene_biotype:protein_coding" Homo_sapiens.GRCh38.cdna.all.fa | awk -F"[> ]" '{print $2;}' > protein_coding.ids
module load GNU/4.4.5
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp Homo_sapiens.GRCh38.cdna.all.fa --output_fasta_fp Homo_sapiens.GRCh38.cdna.ptn.fa --seq_id_fp protein_coding.ids
module load GNU/4.8.3
module load BLAST+/2.3.0
makeblastdb -in Homo_sapiens.GRCh38.cdna.ptn.fa -dbtype nucl -out hg38_cdnaPtn_db

blastn -query novel_I_f2.fa  -db hg38_cdnaPtn_db  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > novel_I_hg38_cdna.outfmt6
blastn -query novel_II_f2.fa  -db hg38_cdnaPtn_db  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > novel_II_hg38_cdna.outfmt6
blastn -query novel_III_f2.fa  -db hg38_cdnaPtn_db  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > novel_III_hg38_cdna.outfmt6
blastn -query intergenic_f2.fa  -db hg38_cdnaPtn_db  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > intergenic_hg38_cdna.outfmt6
blastn -query known_ncRNA_f2.fa  -db hg38_cdnaPtn_db  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > known_lncRNA_hg38_cdna.outfmt6

