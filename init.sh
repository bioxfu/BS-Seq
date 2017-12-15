source activate gmatic
conda env export > doc/environment.yml

if [ ! -d fastq ]; then
	mkdir -p bam clean fastq fastqc/raw fastqc/clean output stat track figures RData
fi
