INPUT=$1
OUTPUT=$2

grep -v '#' $INPUT|grep -v 'log10'|grep -v '^$'|awk '{if($9 < 5)print $1"\t"$2-1"\t"$3"\tMACS_peak_"NR"\t"$7}' > $OUTPUT
