PEAK=peak_1e-2
IP1=IP_ING4.1
IP2=IP_ING5.2
Input1=Input_ING1.1
Input2=Input_ING2.2

cat $PEAK/${IP1}_c_${Input1}_peaks.bed $PEAK/${IP2}_c_${Input2}_peaks.bed|grep -v 'Chr[MC]'|sortBed|mergeBed -i - -c 4,5 -o count,max > $PEAK/peaks_olp_count.tsv

bedtools intersect -a $PEAK/peaks_olp_count.tsv -b bam/${IP1}.bam -c > $PEAK/peaks_olp_count_IP1
bedtools intersect -a $PEAK/peaks_olp_count.tsv -b bam/${IP2}.bam -c > $PEAK/peaks_olp_count_IP2

paste $PEAK/peaks_olp_count_IP1 $PEAK/peaks_olp_count_IP2 |cut -f1-6,12 > $PEAK/peaks_olp_count_coverage

