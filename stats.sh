
known_ncRNA="RNAseqSupTrans.merge.reduced.ncRNA"
P_final_bed <- P_final_bed[!duplicated(P_final_bed), ]
novel_I="inputs/novelAnn/sup/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.sup"
novel_II="inputs/novelAnn/unsup.cons/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.cons"
novel_III="inputs/novelAnn/unsup.uncons.ORF/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.uncons.ORF"
intergenic="inputs/intergenic/unsup.ids"
echo "group" "input_trans" "F1" "F2" "F3" > stats.txt
for gp in known_ncRNA novel_I novel_II novel_III intergenic;do 
 echo "$gp" "$(tail -n+2 ${!gp} | wc -l)" "$(cat F1_$gp | wc -l)" "$(cat F2_$gp | wc -l)" "$(cat ${gp}_P.bed | wc -l)" "$(cat F4_$gp | wc -l)";
done >> stats.txt


echo "group" "input_trans" "af_F1" "af_F2" "af_F3" > stats2.txt
for gp in known_ncRNA novel_I novel_II novel_III intergenic;do 
 echo "$gp" "$(tail -n+2 ${!gp} | wc -l)" "no_file" "$(cat ${gp}_f2.bed | wc -l)" "$(cat ${gp}_f3.bed | wc -l)" "$(cat ${gp}_final.bed | wc -l)";
done >> stats2.txt


echo "rescued_ncRNA:" $(cat lncRNA_rescued.bed | wc -l) >> stats2.txt
echo "total final ncrRNA:" $(cat lncRNA_final.bed | wc -l) >> stats2.txt

