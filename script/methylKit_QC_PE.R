library(methylKit)
library(yaml)

config <- yaml.load_file('config.yaml')
config$samples <- lapply(config$samples, as.character)
context_list <- c('CpG', 'CHH', 'CHG')

for (context in context_list) {
  cat('check', context, '...\n')
  # read the methylation call files
  file_list <- as.list(paste0('output/', config$samples, '_R1_paired_bismark_bt2_pe_sort.deduplicated.', context, '_report.txt.gz'))
  raw <- methRead(file_list, sample.id=config$samples, assembly=config$assembly, treatment=config$treatment, 
                  context=context, pipeline = 'bismarkCytosineReport', header=F)
  
  pdf(paste0('figures/methylKit_QC_', context, '.pdf'))

  # basic stats about the methylation data
  for (i in 1:length(raw)) {
    getMethylationStats(raw[[i]], plot=TRUE)
    getCoverageStats(raw[[i]], plot=TRUE)
  }

  # discards bases that have coverage below 10X 
  # and also discards the bases that have more than 
  # 99.9th percentile of coverage in each sample.
  filtered <- filterByCoverage(raw, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
  
  # merge all samples to one object for base-pair locations that are covered in all samples. 
  meth <- unite(filtered)
  
  # check the correlation between samples
  getCorrelation(meth, plot=TRUE)

  # cluster the samples based on the similarity of their methylation profiles
  clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

  # PCA analysis
  PCASamples(meth)

  dev.off()
}
