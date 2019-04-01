source activate gmatic
#conda env export > doc/environment.yml
#module add macs/1.4.2

if [ ! -d fastqc ]; then
	mkdir fastqc clean bam MAnorm track peak_1e-5 peak_1e-2 table figure RData matrix stat meme-chip
fi

export LD_LIBRARY_PATH=$HOME/miniconda2/envs/gmatic/jre/lib/amd64/server

