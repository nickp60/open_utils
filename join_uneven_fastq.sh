#!/bin/bash
# from bioinformatics magician level 110 Heng Li
# version 0.0.1
# argument 1 is read 1
# argument 2 is read 2


echo 'USAGE: /path/to/contigs.fasta'



awk '{printf substr($0,1,length-2);getline;printf "\t"$0;getline;getline;print "\t"$0}' "$1" | sort -S 8G -T. > read1_temp.txt
awk '{printf substr($0,1,length-2);getline;printf "\t"$0;getline;getline;print "\t"$0}' "$1" | sort -S 8G -T. > read2_temp.txt
join read1_temp.txt read2_temp.txt | awk '{print $1"\n"$2"\n+\n"$3 > "r1.fq";print $1"\n"$4"\n+\n"$5 > "r2.fq"}'

# rm read1_temp.txt -f
# rm read2_temp.txt -f
