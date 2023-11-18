cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGHV.fasta ../IGKL/IGKLV.fasta >IGV.fasta
cat ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGHJ.fasta.V2 ../IGKL/IGKLJ.fasta >IGJ.fasta
sh ../bin/run.sh IGV.fasta IGJ.fasta V310J38 t.txt IG . 1 ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGHD.fasta
