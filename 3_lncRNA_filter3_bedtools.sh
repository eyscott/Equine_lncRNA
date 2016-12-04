##now use bedtools
module load BEDTools/2.24.0  
bedtools intersect -c -a novel_I_f2.bed -b merged_5.bed -sorted -s > novel_I_5.bed
bedtools intersect -c -a novel_II_f2.bed -b merged_5.bed -sorted -s > novel_II_5.bed
bedtools intersect -c -a novel_III_f2.bed -b merged_5.bed -sorted -s > novel_III_5.bed
bedtools intersect -c -a intergenic_f2.bed -b merged_5.bed -sorted -s > intergenic_5.bed


bedtools intersect -c -a novel_I_f2.bed -b merged_3.bed -sorted -s > novel_I_3.bed
bedtools intersect -c -a novel_II_f2.bed -b merged_3.bed -sorted -s > novel_II_3.bed
bedtools intersect -c -a novel_III_f2.bed -b merged_3.bed -sorted -s > novel_III_3.bed
bedtools intersect -c -a intergenic_f2.bed -b merged_3.bed -sorted -s > intergenic_3.bed
