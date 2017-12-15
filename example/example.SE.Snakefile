configfile: "config.yaml"

rule all:
	input:
		expand('clean/{sample}_clean.fastq.gz', sample=config['samples']),
		expand('fastqc/raw/{sample}_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_clean_fastqc.html', sample=config['samples']),
		expand('stat/fastqc_stat.tsv'),
		expand('bam/{sample}_clean_bismark_bt2.deduplicated.bam', sample=config['samples']),
		expand('output/{sample}_clean_bismark_bt2.deduplicated.bedGraph.gz', sample=config['samples']),
		expand('output/{sample}_clean_bismark_bt2.deduplicated.bismark.cov.gz', sample=config['samples']),
		expand('output/{sample}_clean_bismark_bt2.deduplicated.CX_report.txt.CX_report.txt.gz', sample=config['samples']),
		expand('stat/{sample}_clean_bismark_bt2_SE_report.html', sample=config['samples']),
		expand('track/{sample}.tdf', sample=config['samples']),
		expand('output/{sample}_clean_bismark_bt2.deduplicated.CpG_report.txt.gz', sample=config['samples']),
		expand('output/{sample}_clean_bismark_bt2.deduplicated.CHG_report.txt.gz', sample=config['samples']),
		expand('output/{sample}_clean_bismark_bt2.deduplicated.CHH_report.txt.gz', sample=config['samples']),
		expand('figures/methylKit_QC_CpG.pdf'),
		expand('figures/methylKit_QC_CHG.pdf'),
		expand('figures/methylKit_QC_CHH.pdf'),
		expand('RData/methylKit_DMR.RData'),
		
rule fastqc_raw_SE:
	input:
		config['path']+'/{sample}.fastq.gz'
	output:
		'fastqc/raw/{sample}_fastqc.html'
	shell:
		'fastqc -o fastqc/raw {input}'

rule trimmomatic_SE:
	input:
		config['path']+'/{sample}.fastq.gz'
	output:
		'clean/{sample}_clean.fastq.gz'
	params:
		adapter = config['adapter']
	shell:
		'trimmomatic SE -threads 3 -phred33 {input} {output} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule fastqc_clean_SE:
	input:
		'clean/{sample}_clean.fastq.gz'
	output:
		'fastqc/clean/{sample}_clean_fastqc.html'
	shell:
		'fastqc -o fastqc/clean {input}'

rule fastqc_stat_SE:
	input:
		['fastqc/raw/{sample}_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_clean_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule bismark:
	input:
		'clean/{sample}_clean.fastq.gz'
	output:
		'bam/{sample}_clean_bismark_bt2.bam',
		'bam/{sample}_clean_bismark_bt2_SE_report.txt'
	params:
		genome_folder = config['genome_folder'],
		cpu = config['cpu']
	shell:
		'bismark --parallel {params.cpu} --genome {params.genome_folder} -o bam {input}'

rule deduplicate_bismark:
	input:
		'bam/{sample}_clean_bismark_bt2.bam'
	output:
		'bam/{sample}_clean_bismark_bt2.deduplicated.bam',
		'bam/{sample}_clean_bismark_bt2.deduplication_report.txt'
	shell:
		'deduplicate_bismark -s {input} --bam --output_dir bam'

rule bismark_methylation_extractor:
	input:
		'bam/{sample}_clean_bismark_bt2.deduplicated.bam'
	output:
		bedgraph = 'output/{sample}_clean_bismark_bt2.deduplicated.bedGraph.gz',
		cov = 'output/{sample}_clean_bismark_bt2.deduplicated.bismark.cov.gz',
		report = 'output/{sample}_clean_bismark_bt2.deduplicated.CX_report.txt.CX_report.txt.gz',
		splitting = 'output/{sample}_clean_bismark_bt2.deduplicated_splitting_report.txt',
		mbias = 'output/{sample}_clean_bismark_bt2.deduplicated.M-bias.txt',
	params:
		genome_folder = config['genome_folder'],
		cpu = config['cpu']
	shell:
		'bismark_methylation_extractor --comprehensive --CX --gzip --multicore {params.cpu} --bedGraph --cytosine_report --genome_folder {params.genome_folder} -o output {input}'

rule bam2nuc:
	input:
		'bam/{sample}_clean_bismark_bt2.deduplicated.bam'
	output:
		'output/{sample}_clean_bismark_bt2.deduplicated.nucleotide_stats.txt'
	params:
		genome_folder = config['genome_folder']
	shell:
		'bam2nuc --dir output --genome_folder {params.genome_folder} {input}'

rule bismark2report:
	input:
		alignment = 'bam/{sample}_clean_bismark_bt2_SE_report.txt',
		dedup = 'bam/{sample}_clean_bismark_bt2.deduplication_report.txt',
		splitting = 'output/{sample}_clean_bismark_bt2.deduplicated_splitting_report.txt',
		mbias = 'output/{sample}_clean_bismark_bt2.deduplicated.M-bias.txt',
		nucleotide = 'output/{sample}_clean_bismark_bt2.deduplicated.nucleotide_stats.txt'
	output:
		'stat/{sample}_clean_bismark_bt2_SE_report.html'
	shell:
		'bismark2report --dir stat --alignment_report {input.alignment} --dedup_report {input.dedup} --splitting_report {input.splitting} --mbias_report {input.mbias} --nucleotide_report {input.nucleotide}'

rule bedgraph2tdf:
	input:
		bg = 'output/{sample}_clean_bismark_bt2.deduplicated.bedGraph.gz'
	output:
		temp = 'track/{sample}.bedGraph',
		tdf = 'track/{sample}.tdf'
	params:
		IGV = config['IGV']
	shell:
		"zcat {input.bg} > {output.temp}; igvtools toTDF {output.temp} {output.tdf} {params.IGV}"

rule context_type:
	input:
		'output/{sample}_clean_bismark_bt2.deduplicated.CX_report.txt.CX_report.txt.gz'
	output:
		CpG = 'output/{sample}_clean_bismark_bt2.deduplicated.CpG_report.txt.gz',
		CHG = 'output/{sample}_clean_bismark_bt2.deduplicated.CHG_report.txt.gz',
		CHH = 'output/{sample}_clean_bismark_bt2.deduplicated.CHH_report.txt.gz'
	shell:
		"zcat {input}|grep -v 'CH[HG]'|gzip -c > {output.CpG}; zcat {input}|grep 'CHG'|gzip -c > {output.CHG}; zcat {input}|grep 'CHH'|gzip -c > {output.CHH}"

rule methylKit_QC:
	input:
		['output/{sample}_clean_bismark_bt2.deduplicated.CpG_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_clean_bismark_bt2.deduplicated.CHG_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_clean_bismark_bt2.deduplicated.CHH_report.txt.gz'.format(sample=x) for x in config['samples']]
	output:
		'figures/methylKit_QC_CpG.pdf',
		'figures/methylKit_QC_CHG.pdf',
		'figures/methylKit_QC_CHH.pdf'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/methylKit_QC.R'

rule methylKit_DMR:
	input:
		['output/{sample}_clean_bismark_bt2.deduplicated.CpG_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_clean_bismark_bt2.deduplicated.CHG_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_clean_bismark_bt2.deduplicated.CHH_report.txt.gz'.format(sample=x) for x in config['samples']]
	output:
		'RData/methylKit_DMR.RData'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/methylKit_DMR.R'

