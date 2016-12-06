##now use bedtools
module load GNU/4.4.5 
module load BEDTools/2.24.0  
bedtools intersect -c -a novel_I_f3.bed -b refined_nolncRNA_5.bed -sorted -s -f 0.25 > novel_I_5.bed
bedtools intersect -c -a novel_II_f3.bed -b refined_nolncRNA_5.bed -sorted -s -f 0.25 > novel_II_5.bed
bedtools intersect -c -a novel_III_f3.bed -b refined_nolncRNA_5.bed -sorted -s -f 0.25 > novel_III_5.bed
bedtools intersect -c -a intergenic_f3.bed -b refined_nolncRNA_5.bed -sorted -s -f 0.25 > intergenic_5.bed
bedtools intersect -c -a known_lncRNA_f3.bed -b refined_nolncRNA_5.bed -sorted -s -f 0.25 > known_5.bed

bedtools intersect -c -a novel_I_f3.bed -b refined_nolncRNA_3.bed -sorted -s -f 0.25 > novel_I_3.bed
bedtools intersect -c -a novel_II_f3.bed -b refined_nolncRNA_3.bed -sorted -s -f 0.25 > novel_II_3.bed
bedtools intersect -c -a novel_III_f3.bed -b refined_nolncRNA_3.bed -sorted -s -f 0.25 > novel_III_3.bed
bedtools intersect -c -a intergenic_f3.bed -b refined_nolncRNA_3.bed -sorted -s -f 0.25 > intergenic_3.bed
bedtools intersect -c -a known_lncRNA_f3.bed -b refined_nolncRNA_3.bed -sorted -s -f 0.25 > known_3.bed
