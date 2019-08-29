IP=$1
INPUT=$2

grep 'Chr1' $IP > $IP.Chr1
grep 'Chr2' $IP > $IP.Chr2
grep 'Chr3' $IP > $IP.Chr3
grep 'Chr4' $IP > $IP.Chr4
grep 'Chr5' $IP > $IP.Chr5
grep 'ChrC' $IP > $IP.ChrC
grep 'ChrM' $IP > $IP.ChrM
 
grep 'Chr1' $INPUT > $INPUT.Chr1
grep 'Chr2' $INPUT > $INPUT.Chr2
grep 'Chr3' $INPUT > $INPUT.Chr3
grep 'Chr4' $INPUT > $INPUT.Chr4
grep 'Chr5' $INPUT > $INPUT.Chr5
grep 'ChrC' $INPUT > $INPUT.ChrC
grep 'ChrM' $INPUT > $INPUT.ChrM
 
bedtools intersect -a $IP.Chr1 -b $INPUT.Chr1 -v >  $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr2 -b $INPUT.Chr2 -v >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr3 -b $INPUT.Chr3 -v >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr4 -b $INPUT.Chr4 -v >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr5 -b $INPUT.Chr5 -v >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.ChrC -b $INPUT.ChrC -v >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.ChrM -b $INPUT.ChrM -v >> $IP.subtract_Input.bedgraph
 
bedtools intersect -a $IP.Chr1 -b $INPUT.Chr1 -wa -wb|cut -f1-4,8|groupBy -g 1,2,3,4 -c 5 -o mean|awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr2 -b $INPUT.Chr2 -wa -wb|cut -f1-4,8|groupBy -g 1,2,3,4 -c 5 -o mean|awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr3 -b $INPUT.Chr3 -wa -wb|cut -f1-4,8|groupBy -g 1,2,3,4 -c 5 -o mean|awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr4 -b $INPUT.Chr4 -wa -wb|cut -f1-4,8|groupBy -g 1,2,3,4 -c 5 -o mean|awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.Chr5 -b $INPUT.Chr5 -wa -wb|cut -f1-4,8|groupBy -g 1,2,3,4 -c 5 -o mean|awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.ChrC -b $INPUT.ChrC -wa -wb|cut -f1-4,8|groupBy -g 1,2,3,4 -c 5 -o mean|awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' >> $IP.subtract_Input.bedgraph
bedtools intersect -a $IP.ChrM -b $INPUT.ChrM -wa -wb|cut -f1-4,8|groupBy -g 1,2,3,4 -c 5 -o mean|awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' >> $IP.subtract_Input.bedgraph
 
sortBed -i $IP.subtract_Input.bedgraph > $IP.subtract_Input.bedgraph.tmp; mv $IP.subtract_Input.bedgraph.tmp $IP.subtract_Input.bedgraph
rm $IP.Chr* $INPUT.Chr*

igvtools toTDF $IP.subtract_Input.bedgraph $IP.subtract_Input.tdf $HOME/igv/genomes/tair10.genome
