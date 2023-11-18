 cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRAV.fasta ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRBV.fasta >TRABV.fasta
cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRAJ.fasta.V2 ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRBJ.fasta.V2 >TRABJ.fasta
sh ../bin/run.sh TRABV.fasta TRABJ.fasta V310J38 t.txt TRAB . 1 ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRBD.fasta
