cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRDV.fasta ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRGV.fasta >TRDGV.fasta
cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRDJ.fasta.V2 ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRGJ.fasta.V2 >>TRDGJ.fasta
sh ../bin/run.sh TRDGV.fasta TRDGJ.fasta V310J38 t.txt TRDG . 1 ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/TR_format/TRDD.fasta
