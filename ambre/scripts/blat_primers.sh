
#~/software/blat -t=dna -q=dna -stepSize=2 -tileSize=8 -repMatch=1048576 -minScore=15 -noHead ~/data/hg/hg19/full.fa gbm-cDNA01.a1.sln00.fa gbm-cDNA01.a1.sln00.blat.hg19.fa
sort -nr gbm-cDNA03.a3.sln00.blat.hg19.out | awk '{if($8<2 && ($14=="chr4" || $14=="chr7")){print $14 "\t" $16+1 "\t" $17+1 "\t" $10;}}' | python quick_unique.py 4 >gbm-cDNA03.a3.sln00.blat.hg19.bed
