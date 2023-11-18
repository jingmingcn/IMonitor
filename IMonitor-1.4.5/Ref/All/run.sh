cat ../TR/TRV.fasta ../IG/IGV.fasta >AllV.fasta
cat ../TR/TRJ.fasta ../IG/IGJ.fasta >AllJ.fasta
cat ../TR/TRD.fasta ../../../IMGT_V-QUEST_reference_directory_20210531/Homo_sapiens/IG_format/IGHD.fasta >AllD.fasta
sh ../bin/run.sh AllV.fasta AllJ.fasta V310J38 t.txt All . 1 AllD.fasta
