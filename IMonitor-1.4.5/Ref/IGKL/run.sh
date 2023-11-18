cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGKV.fasta ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGLV.fasta >IGKLV.fasta
cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGKJ.fasta.V2 ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGLJ.fasta.V2 >IGKLJ.fasta
sh ../bin/run.sh IGKLV.fasta IGKLJ.fasta V310J38 t.txt IGKL . 1
