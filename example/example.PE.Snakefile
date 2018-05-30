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
		expand('bam/{sample}_bismark_sort.deduplicated.bam', sample=config['samples']),
		expand('output/{sample}_bismark_sort.deduplicated.CX_report.txt.gz', sample=config['samples']),
		expand('stat/{sample}_bismark_report.html', sample=config['samples']),
		expand('output/{sample}_bismark_sort.deduplicated.{context}_report.txt.gz', sample=config['samples'], context=['CpG', 'CHG', 'CHH']),
		expand('track/{sample}_{context}_min{min}.tdf', sample=config['samples'], context=['CpG', 'CHG', 'CHH'], min=config['min_coverage']),
		['tables/methyC_{context}_counts.tsv'.format(context=x) for x in ['CpG', 'CHG', 'CHH']],
		['tables/methyC_{context}_counts_filter.tsv'.format(context=x) for x in ['CpG', 'CHG', 'CHH']],
		['figures/methyC_{context}_counts_PCA.pdf'.format(context=x) for x in ['CpG', 'CHG', 'CHH']],

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

rule rename_bismark_output_bam:
	input:
		bam = 'bam/{sample}_R1_paired_bismark_bt2_pe.bam',
		report = 'bam/{sample}_R1_paired_bismark_bt2_PE_report.txt'
	output:
		bam = 'bam/{sample}_bismark.bam',
		report = 'bam/{sample}_bismark_report.txt'
	shell:
		'mv {input.bam} {output.bam}; mv {input.report} {output.report}'

rule sort_bam_by_name:
	input:
		'bam/{sample}_bismark.bam',
	output:
		'bam/{sample}_bismark_sort.bam',
	shell:
		'samtools sort -n {input} -o {output}'

rule deduplicate_bismark:
	input:
		'bam/{sample}_bismark_sort.bam',
	output:
		'bam/{sample}_bismark_sort.deduplicated.bam',
		'bam/{sample}_bismark_sort.deduplication_report.txt'
	shell:
		'deduplicate_bismark -p {input} --bam --output_dir bam'

rule bam2nuc:
	input:
		'bam/{sample}_bismark_sort.deduplicated.bam'
	output:
		'output/{sample}_bismark_sort.deduplicated.nucleotide_stats.txt'
	params:
		genome_folder = config['genome_folder']
	shell:
		'bam2nuc --dir output --genome_folder {params.genome_folder} {input}'

rule bismark_methylation_extractor:
	input:
		'bam/{sample}_bismark_sort.deduplicated.bam'
	output:
		report = 'output/{sample}_bismark_sort.deduplicated.CX_report.txt.CX_report.txt.gz',
		splitting = 'output/{sample}_bismark_sort.deduplicated_splitting_report.txt',
		mbias = 'output/{sample}_bismark_sort.deduplicated.M-bias.txt',
	params:
		genome_folder = config['genome_folder'],
		cpu = config['cpu']
	shell:
		'bismark_methylation_extractor --comprehensive --CX --gzip --multicore {params.cpu} --cytosine_report --genome_folder {params.genome_folder} -o output {input}'

rule rename_bismark_output_report:
	input:
		report = 'output/{sample}_bismark_sort.deduplicated.CX_report.txt.CX_report.txt.gz',
	output:
		report = 'output/{sample}_bismark_sort.deduplicated.CX_report.txt.gz',
	shell:
		'mv {input.report} {output.report}'

rule bismark2report:
	input:
		alignment = 'bam/{sample}_bismark_report.txt',
		dedup = 'bam/{sample}_bismark_sort.deduplication_report.txt',
		splitting = 'output/{sample}_bismark_sort.deduplicated_splitting_report.txt',
		mbias = 'output/{sample}_bismark_sort.deduplicated.M-bias.txt',
		nucleotide = 'output/{sample}_bismark_sort.deduplicated.nucleotide_stats.txt'
	output:
		'stat/{sample}_bismark_report.html'
	shell:
		'bismark2report --dir stat --alignment_report {input.alignment} --dedup_report {input.dedup} --splitting_report {input.splitting} --mbias_report {input.mbias} --nucleotide_report {input.nucleotide}'

rule context_type_report:
	input:
		'output/{sample}_bismark_sort.deduplicated.CX_report.txt.gz'
	output:
		CpG = 'output/{sample}_bismark_sort.deduplicated.CpG_report.txt.gz',
		CHG = 'output/{sample}_bismark_sort.deduplicated.CHG_report.txt.gz',
		CHH = 'output/{sample}_bismark_sort.deduplicated.CHH_report.txt.gz'
	shell:
		"zcat {input}|grep -v 'CH[HG]'|cut -f1-5|sort -k1,1 -k2,2n|gzip -c > {output.CpG}; zcat {input}|grep 'CHG'|cut -f1-5|sort -k1,1 -k2,2n|gzip -c > {output.CHG}; zcat {input}|grep 'CHH'|cut -f1-5|sort -k1,1 -k2,2n|gzip -c > {output.CHH}"

rule report2bedgraph:
	input: 
		'output/{sample}_bismark_sort.deduplicated.{context}_report.txt.gz',
	output:
		bedgraph = 'track/{sample}_{context}_min{min}.bedGraph',
		tdf = 'track/{sample}_{context}_min{min}.tdf'
	params:
		IGV = config['IGV'],
		min_coverage = config['min_coverage']
	shell:
		'script/methy_bedgraph.py {params.min_coverage} {input} > {output.bedgraph}; igvtools toTDF {output.bedgraph} {output.tdf} {params.IGV}'

rule merge_CX_report_table:
	input:
		['output/{sample}_bismark_sort.deduplicated.CpG_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_bismark_sort.deduplicated.CHH_report.txt.gz'.format(sample=x) for x in config['samples']],
		['output/{sample}_bismark_sort.deduplicated.CHG_report.txt.gz'.format(sample=x) for x in config['samples']]
	output:
		'tables/methyC_CpG_counts.tsv',
		'tables/methyC_CHG_counts.tsv',
		'tables/methyC_CHH_counts.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/merge_CX_report_table.R'

rule min_coverage_filter:
	input:
		'tables/methyC_{context}_counts.tsv'
	output:
		'tables/methyC_{context}_counts_filter.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/min_coverage_filter.R {input} {output}'

rule methy_PCA:
	input:
		'tables/methyC_{context}_counts.tsv'
	output:
		'figures/methyC_{context}_counts_PCA.pdf'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/methy_PCA.R {input} {output}'
