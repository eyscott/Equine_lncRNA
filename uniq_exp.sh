cp inputs/backmapping_stats/allTissues_isoformTPM allTrans_allTissues_isoformTPM
tail -n+2 inputs/backmapping_stats/intergenic_allTissues_isoformTPM >> allTrans_allTissues_isoformTPM

head -n+1 allTrans_allTissues_isoformTPM > lncRNA_allTissues_isoformTPM
awk -F"\t" '{print $4}' lncRNA_final.bed | grep -Fwf - allTrans_allTissues_isoformTPM  >> lncRNA_allTissues_isoformTPM

targets=()
for target in BrainStem Cerebellum Embryo.ICM Embryo.TE Muscle Retina Skin SpinalCord;do targets+=($target);done

##print no of gene/isoform expressed, expressed uniqely, not expressed uniquely
mkdir -p uniqExp_lncRNA
cutoff=0.1            ## run with cutoff=0or5
i=0
n=${#targets[@]}
echo "## no of gene/isoform expressed, expressed uniqely, not expressed uniquely" > tissueSpecificSummary_cutoff.$cutoff
while [ $i -lt $n ];do
 echo ${targets[$i]}
 #cat allTissues_geneTPM | awk -v x=$((i+2)) -v c=$cutoff '$x>c' | wc -l
 cat lncRNA_allTissues_isoformTPM | awk -v x=$((i+2)) -v c=$cutoff '$x>c' | wc -l

 #cat allTissues_geneTPM > tempGene;
 cat lncRNA_allTissues_isoformTPM > tempIsoform;
 for x in `seq 2 $((n+1))`;do
  #if [ $x -eq $((i+2)) ];then awk -v x=$x -v c=$cutoff '$x>c' tempGene > tempGene2; else awk -v x=$x -v c=$cutoff '$x<=c' tempGene > tempGene2;fi
  if [ $x -eq $((i+2)) ];then awk -v x=$x -v c=$cutoff '$x>c' tempIsoform > tempIsoform2; else awk -v x=$x -v c=$cutoff '$x<=c' tempIsoform > tempIsoform2;fi
  #mv tempGene2 tempGene;
  mv tempIsoform2 tempIsoform;done
 #cat tempGene | wc -l;
 cat tempIsoform | wc -l;
 #mv tempGene uniqExp_lncRNA/${targets[$i]}.gene.expressed_uniqely_cutoff.$cutoff;
 mv tempIsoform uniqExp_lncRNA/${targets[$i]}.isoform.expressed_uniqely_cutoff.$cutoff;

 #cat allTissues_geneTPM > tempGene;
 cat lncRNA_allTissues_isoformTPM > tempIsoform;
 for x in `seq 2 $((n+1))`;do
  #if [ $x -eq $((i+2)) ];then awk -v x=$x -v c=$cutoff '$x<=c' tempGene > tempGene2; else awk -v x=$x -v c=$cutoff '$x>c' tempGene > tempGene2;fi
  if [ $x -eq $((i+2)) ];then awk -v x=$x -v c=$cutoff '$x<=c' tempIsoform > tempIsoform2; else awk -v x=$x -v c=$cutoff '$x>c' tempIsoform > tempIsoform2;fi
  #mv tempGene2 tempGene;
  mv tempIsoform2 tempIsoform;done
 #cat tempGene | wc -l;
 cat tempIsoform | wc -l;
 #mv tempGene uniqExp_lncRNA/${targets[$i]}.gene.notExpressed_uniqely_cutoff.$cutoff;
 mv tempIsoform uniqExp_lncRNA/${targets[$i]}.isoform.notExpressed_uniqely_cutoff.$cutoff;
 ((i+=1))
done >> tissueSpecificSummary_cutoff.$cutoff

