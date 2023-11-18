cat ../TRAB/TRABV.fasta ../TRDG/TRDGV.fasta >TRV.fasta
cat ../TRAB/TRABJ.fasta ../TRDG/TRDGJ.fasta >TRJ.fasta
cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRBD.fasta ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRDD.fasta >TRD.fasta
sh ../bin/run.sh TRV.fasta TRJ.fasta V310J38 t.txt TR . 1 TRD.fasta
