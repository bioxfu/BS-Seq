source activate gmatic
conda env export > doc/environment.yml
module add bcftools/0.1.19
module add samtools/0.1.19

if [ ! -d fastq ]; then
	mkdir -p bam clean fastq fastqc/raw fastqc/clean output stat track figures RData
fi
