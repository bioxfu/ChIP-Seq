ls ./fastq/*/*/*.gz |./script/rush -k 'cat {} >> fastq/{/%}_R{%@.+[_.]R*([12]).fa*s*t*q.gz}.fastq.gz' --dry-run | tee combine.sh
echo 'combine multiple runs for each sample...'
bash combine.sh
echo 'done'
