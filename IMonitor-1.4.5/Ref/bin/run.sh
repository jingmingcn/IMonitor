if [ $# -ne 8 -a $# -ne 7 ]
then
	echo "USAGE:  sh $0 <V_sort.fa> <J_sort.fa> <cdr3region> <primer.txt> <gene_type> <out_dir> <flag> [<D_sort.fa>]"
	echo -e "\tIf create the normal reference, then <flag> = 1 && <primer.txt> is invail and use any file input is OK\n"
	echo -e "\t<V_sort.fa>/<J_sort.fa>/<D_sort.fa>: V/D/J germline FASTA sequences, with sorted\n";
	echo -e "\t<cdr3region>: the position of CDR3 start in V and CDR3 end in J. e.g. V310J25, means CDR3 start from the 310the position of V, end at the 25th position of J\n";
	echo -e "\t<gene_type>: gene name, e.g. IGH,TRB,TRA,IGKL\n";
	exit 
else
	V_fa=$1
	J_fa=$2
	cdr3region=$3
	primer=$4
	gene_type=$5
	refdir=$6
	flag=$7
	D_fa=$8
fi

#bin="/hwfssz1/ST_HEALTH/Immune_And_Health_Lab/Public_Software/Immune_repertoire/Common_used/IMonitor-1.3.1/Ref/bin/"
bin="/data/Public_tools/Pipeline/pipeline_hwfssz1/IMonitor/IMonitor-1.4.1/Ref/bin/"

# step1: find primer: Multiplex and 5'RACE use different program ------------------------
if [ $flag -eq 1 ]
then
	perl $bin/5race_V_seq_format.pl $V_fa 1 $refdir/${gene_type}V_process_primer_ref_sort.fa 
	perl $bin/5race_V_seq_format.pl $J_fa 1 $refdir/${gene_type}J_process_primer_ref_sort.fa
else
	perl $bin/Use_primer_find_pcr_seq.pl $primer $V_fa V $refdir/${gene_type}V_process_primer_ref_sort.fa
	perl $bin/Use_primer_find_pcr_seq.pl $primer $J_fa J $refdir/${gene_type}J_process_primer_ref_sort.fa
fi


# step2: combine sequences: judge if exists D gene --------------------------
if [ -e "$D_fa" ]
then
	perl $bin/Ref_seq_integration.pl $refdir/${gene_type}V_process_primer_ref_sort.fa  $D_fa $refdir/${gene_type}J_process_primer_ref_sort.fa $cdr3region $gene_type $refdir >$refdir/${gene_type}.primer.backup
	cat $refdir/${gene_type}.V.alignment $refdir/${gene_type}.D.alignment $refdir/${gene_type}.J.alignment >$refdir/${gene_type}.alignment
	rm $refdir/${gene_type}.V.alignment $refdir/${gene_type}.D.alignment $refdir/${gene_type}.J.alignment
else
	perl $bin/Ref_seq_integration.pl $refdir/${gene_type}V_process_primer_ref_sort.fa $refdir/${gene_type}J_process_primer_ref_sort.fa $cdr3region $gene_type $refdir >$refdir/${gene_type}.primer.backup
	 cat $refdir/${gene_type}.V.alignment $refdir/${gene_type}.J.alignment >$refdir/${gene_type}.alignment
	 rm $refdir/${gene_type}.V.alignment $refdir/${gene_type}.J.alignment
fi

rm $refdir/${gene_type}V_process_primer_ref_sort.fa $refdir/${gene_type}J_process_primer_ref_sort.fa 
rm $refdir/${gene_type}.converted.SingleLine



# step3: blast index ---------------------------
if [ -e "$D_fa" ]
then
	$bin/formatdb -i $refdir/$gene_type.D.converted.fa -p F -o T -e T -n $refdir/$gene_type.D.ref_index
fi

$bin/formatdb -i $refdir/$gene_type.V.converted.fa -p F -o T -e T -n $refdir/$gene_type.V.ref_index
$bin/formatdb -i $refdir/$gene_type.J.converted.fa -p F -o T -e T -n $refdir/$gene_type.J.ref_index

