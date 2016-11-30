#step 1: download Pfam-A db from:  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

#step 2: get the fasta format of the lncRNA from the .bed file:
#using bedtools to get the fasta
#getFastaFromBed [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>

##main script for all novel categories passing the filters so far
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab2/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
tar xvzf chromFa.tar.gz
cat *.fa > equCab2.0_genome.fa
module load BEDTools/2.24.0  
bedtools getfasta -s -fi equCab2.0_genome.fa -bed novel_I_f3_new.bed -fo novel_I_f3.fa
bedtools getfasta -s -fi equCab2.0_genome.fa -bed novel_II_f3_new.bed -fo novel_II_f3.fa
bedtools getfasta -s -fi equCab2.0_genome.fa -bed novel_III_f3_new.bed -fo novel_III_f3.fa
bedtools getfasta -s -fi equCab2.0_genome.fa -bed intergenic_f3_new.bed -fo intergenic_f3.fa

#step 3: get the protein coding/ORF sequences of your data using Transdecoder (this takes quite some time):
#main
module load TransDecoder/3.0.0
TransDecoder.LongOrfs -t novel_I_f3.fa -m 20
TransDecoder.LongOrfs -t novel_II_f3.fa -m 20
TransDecoder.LongOrfs -t novel_III_f3.fa -m 20
TransDecoder.LongOrfs -t intergenic_f3.fa -m 20 
#with -m 20, meaning ORF minimum of 20 (default=100) because we are dealing with lncRNA

##step 4: Predict likely coding regions with HMMER: hmmsearch, using ORFs predicted above:
module load HMMER/3.1b2
#main
hmmsearch --cpu 8 --domtblout novel_I_pfam.domtblout Pfam-A.hmm novel_I_f3.fa.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu 8 --domtblout novel_II_pfam.domtblout Pfam-A.hmm novel_II_f3.fa.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu 8 --domtblout novel_III_pfam.domtblout Pfam-A.hmm novel_III_f3.fa.transdecoder_dir/longest_orfs.pep
hmmsearch --cpu 8 --domtblout intergenic_pfam.domtblout Pfam-A.hmm intergenic_f3.fa.transdecoder_dir/longest_orfs.pep

##step 5: retain those significant hits to the Pfam-A db using Transdecoder:
#main
#TransDecoder.Predict -t novel_I.fa --retain_pfam_hits novel_I_pfam.domtblout
#TransDecoder.Predict -t novel_II.fa --retain_pfam_hits novel_II_pfam.domtblout
#TransDecoder.Predict -t novel_III.fa --retain_pfam_hits novel_III_pfam.domtblout
#TransDecoder.Predict -t intergenic.fa --retain_pfam_hits intergenic_pfam.domtblout

###step 6 running blastp
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
module load GNU/4.8.3
module load BLAST+/2.3.0
blastp -query novel_I_f3.fa.transdecoder_dir/longest_orfs.pep  -db uniref90.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > novel_I_blastp.outfmt6
blastp -query novel_II_f3.fa.transdecoder_dir/longest_orfs.pep  -db uniref90.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > novel_II_blastp.outfmt6
blastp -query novel_III_f3.fa.transdecoder_dir/longest_orfs.pep  -db uniref90.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > novel_III_blastp.outfmt6
blastp -query intergenic_f3.fa.transdecoder_dir/longest_orfs.pep  -db uniref90.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > intergenic_blastp.outfmt6
