zcat clean/IP_5DAP_lec1l1l.3_R1_paired.fastq.gz|grep '^[ATCG]'|cut -c1-50|head -50|awk '{print ">"NR"\n"$0}' > test.fa
zcat clean/Input_5DAP_WT.1_R1_paired.fastq.gz|grep '^[ATCG]'|cut -c1-50|head -50|awk '{print ">"NR"\n"$0}' > test2.fa

./web_blast.pl megablast nt test.fa > test.out
./web_blast.pl megablast nt test2.fa > test2.out

grep '^>' test.out |sed 's/>//'|sed 's/ /\t/'|sort -k2 > test.out.tab
grep '^>' test2.out |sed 's/>//'|sed 's/ /\t/'|sort -k2 > test2.out.tab

