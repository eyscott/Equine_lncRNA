module load GNU/4.4.5
module load BEDTools/2.24.0
#getting the sequences not associating with the above
tr '\r' '\n' < equCab2.chrom.sizes.txt > equCab2.chrom.sizes
bedtools complement -i equine_transcription.bed -g equCab2.chrom.sizes > intergenic.genome


bedtools random -g intergenic.genome -l 1000 -n 5000 > random_EC.bed

bedtools getfasta -s -fi equCab2.0_genome.fa -bed random_EC.bed -name -split -fo random.fa

module load GNU/4.8.3
module load BLAST+/2.3.0

blastn -query random.fa  -db hg38_genome_db  -max_target_seqs 1 -outfmt "6 qseqid sseqid pident qcovs evalue" -evalue 1e-5 -num_threads 10 > random_genome.outfmt6
