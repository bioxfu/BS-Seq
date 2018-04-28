configfile: "config.yaml"

rule all:
	input:
		expand('clean/{sample}_R1_paired.fastq.gz', sample=config['samples']),
		expand('clean/{sample}_R2_paired.fastq.gz', sample=config['samples']),
		expand('fastqc/raw/{sample}_R1_fastqc.html', sample=config['samples']),
		expand('fastqc/raw/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R1_paired_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R2_paired_fastqc.html', sample=config['samples']),
		expand('stat/fastqc_stat.tsv'),
		expand('bam/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.bam', sample=config['samples']),
		expand('stat/{sample}_R1_paired_bismark_bt2_PE_report.html', sample=config['samples']),
		expand('output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CpG_report.txt.gz', sample=config['samples']),
		expand('output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHG_report.txt.gz', sample=config['samples']),
		expand('output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHH_report.txt.gz', sample=config['samples']),
		expand('track/{sample}.bedGraph', sample=config['samples']),
		expand('track/{sample}.tdf', sample=config['samples']),
		expand('figures/methylKit_PCA.pdf'),
#		expand('RData/methylKit_DMR.RData'),
#		expand('tables/methylKit_DMR_anno.tsv'),
		expand('tables/DMC_fisher_test_anno.tsv'),

rule fastqc_raw_PE:
	input:
		config['path']+'/{sample}_R1.fastq.gz',
		config['path']+'/{sample}_R2.fastq.gz'
	output:
		'fastqc/raw/{sample}_R1_fastqc.html',
		'fastqc/raw/{sample}_R2_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/raw {input}'

rule trimmomatic_PE:
	input:
		r1 = config['path']+'/{sample}_R1.fastq.gz',
		r2 = config['path']+'/{sample}_R2.fastq.gz'
	output:
		r1_paired = 'clean/{sample}_R1_paired.fastq.gz',
		r2_paired = 'clean/{sample}_R2_paired.fastq.gz',
		r1_unpaired = 'clean/{sample}_R1_unpaired.fastq.gz',
		r2_unpaired = 'clean/{sample}_R2_unpaired.fastq.gz'
	params:
		adapter = config['adapter']
	shell:
		'trimmomatic PE -threads 3 -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule fastqc_clean_PE:
	input:
		'clean/{sample}_R1_paired.fastq.gz',
		'clean/{sample}_R2_paired.fastq.gz'
	output:
		'fastqc/clean/{sample}_R1_paired_fastqc.html',
		'fastqc/clean/{sample}_R2_paired_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/clean {input}'

rule fastqc_stat_PE:
	input:
		['fastqc/raw/{sample}_R1_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/raw/{sample}_R2_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R1_paired_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R2_paired_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule bismark_PE:
	input:
		r1 = 'clean/{sample}_R1_paired.fastq.gz',
		r2 = 'clean/{sample}_R2_paired.fastq.gz'
	output:
		'bam/{sample}_R1_paired_bismark_bt2_pe.bam',
		'bam/{sample}_R1_paired_bismark_bt2_PE_report.txt'
	params:
		genome_folder = config['genome_folder'],
		cpu = config['cpu']
	shell:
		'bismark --parallel {params.cpu} --genome {params.genome_folder} -o bam -1 {input.r1} -2 {input.r2}'

rule sort_bam_by_name:
	input:
		'bam/{sample}_R1_paired_bismark_bt2_pe.bam',
	output:
		'bam/{sample}_R1_paired_bismark_bt2_pe_sort.bam',
	shell:
		'samtools sort -n {input} -o {output}'

rule deduplicate_bismark:
	input:
		'bam/{sample}_R1_paired_bismark_bt2_pe_sort.bam',
	output:
		'bam/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.bam',
		'bam/{sample}_R1_paired_bismark_bt2_pe_sort.deduplication_report.txt'
	shell:
		'deduplicate_bismark -p {input} --bam --output_dir bam'

rule bismark_methylation_extractor:
	input:
		'bam/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.bam'
	output:
#		bedgraph = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.bedGraph.gz',
#		cov = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.bismark.cov.gz',
		report = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CX_report.txt.CX_report.txt.gz',
		splitting = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated_splitting_report.txt',
		mbias = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.M-bias.txt',
	params:
		genome_folder = config['genome_folder'],
		cpu = config['cpu']
	shell:
		#'bismark_methylation_extractor --comprehensive --CX --gzip --multicore {params.cpu} --bedGraph --cytosine_report --genome_folder {params.genome_folder} -o output {input}'
		'bismark_methylation_extractor --comprehensive --CX --gzip --multicore {params.cpu} --cytosine_report --genome_folder {params.genome_folder} -o output {input}'

rule bam2nuc:
	input:
		'bam/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.bam'
	output:
		'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.nucleotide_stats.txt'
	params:
		genome_folder = config['genome_folder']
	shell:
		'bam2nuc --dir output --genome_folder {params.genome_folder} {input}'

rule bismark2report:
	input:
		alignment = 'bam/{sample}_R1_paired_bismark_bt2_PE_report.txt',
		dedup = 'bam/{sample}_R1_paired_bismark_bt2_pe_sort.deduplication_report.txt',
		splitting = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated_splitting_report.txt',
		mbias = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.M-bias.txt',
		nucleotide = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.nucleotide_stats.txt'
	output:
		'stat/{sample}_R1_paired_bismark_bt2_PE_report.html'
	shell:
		'bismark2report --dir stat --alignment_report {input.alignment} --dedup_report {input.dedup} --splitting_report {input.splitting} --mbias_report {input.mbias} --nucleotide_report {input.nucleotide}'

rule context_type:
	input:
		'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CX_report.txt.CX_report.txt.gz'
	output:
		CpG = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CpG_report.txt.gz',
		CHG = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHG_report.txt.gz',
		CHH = 'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHH_report.txt.gz'
	shell:
		"zcat {input}|grep -v 'CH[HG]'|gzip -c > {output.CpG}; zcat {input}|grep 'CHG'|gzip -c > {output.CHG}; zcat {input}|grep 'CHH'|gzip -c > {output.CHH}"

rule report2bedgraph:
	input: 
		'output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CX_report.txt.CX_report.txt.gz'
	output:
		bedgraph = 'track/{sample}.bedGraph',
		tdf = 'track/{sample}.tdf'
	params:
		IGV = config['IGV']
	shell:
		'script/methy_bedgraph.py {input}|sort -k1,1 -k2,2n > {output.bedgraph}; igvtools toTDF {output.bedgraph} {output.tdf} {params.IGV}'

rule methylKit_PCA_PE:
	input:
		['output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CpG_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHG_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHH_report.txt.gz'.format(sample=x) for x in config['samples']]
	output:
		'figures/methylKit_PCA.pdf',
		'tables/all_CpG_counts.tsv',
		'tables/all_CHH_counts.tsv',
		'tables/all_CHG_counts.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/methylKit_PCA_PE.R'

#rule methylKit_DMR:
#	input:
#		['output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CpG_report.txt.gz'.format(sample=x) for x in config['samples']],
#		['output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHG_report.txt.gz'.format(sample=x) for x in config['samples']],
#		['output/{sample}_R1_paired_bismark_bt2_pe_sort.deduplicated.CHH_report.txt.gz'.format(sample=x) for x in config['samples']]
#	output:
#		'RData/methylKit_DMR.RData'
#	params:
#		Rscript = config['Rscript_path']
#	shell:
#		'{params.Rscript} script/methylKit_DMR_PE.R'

#rule methylKit_DMR_anno:
#	input:
#		'RData/methylKit_DMR.RData'
#	output:
#		'tables/methylKit_DMR_anno.tsv'	
#	params:
#		Rscript = config['Rscript_path']
#	shell:
#		'{params.Rscript} script/methylKit_anno.R'

rule DMC_anno:
	input:
		'tables/all_CpG_counts.tsv',
		'tables/all_CHH_counts.tsv',
		'tables/all_CHG_counts.tsv'
	output:
		'tables/DMC_fisher_test_anno.tsv'	
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/fisher_test.R'
