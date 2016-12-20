module load GNU/4.4.5
module load BEDTools/2.24.0
#getting the sequences not associating with the above
cp equine_transcription.bed equine_transcription_safe.bed
while read chr size;do
 awk -F"\t" -v chr=$chr -v size=$size '{if($1==chr && $3>size)next;print;}' equine_transcription_safe.bed > temp.bed
 cp temp2.bed temp.bed
done < equCab2.chrom.sizes
diff equine_transcription.bed equine_transcription_safe.bed

bedtools complement -i equine_transcription_safe.bed -g equCab2.chrom.sizes > intergenic.genome

bedtools random -g intergenic.genome -l 1000 -n 5000 > random_EC.bed

bedtools getfasta -s -fi equCab2.0_genome.fa -bed random_EC.bed -name -split -fo random.fa
