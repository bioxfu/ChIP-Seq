## run ##
nohup snakemake -j 30 -rp --latency-wait 60 >> nohup.log 2>&1 &
