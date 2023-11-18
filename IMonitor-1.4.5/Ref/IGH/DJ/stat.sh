le D_ref.raw |awk 'NR>1'|awk '{print ">|"$1"||F\n"$2;}'|perl /data/Public_tools/Pipeline/bin/chang_fa_1_line.pl - >IGH_D_ref.fa
