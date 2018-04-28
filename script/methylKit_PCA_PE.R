library(methylKit)
library(yaml)

config <- yaml.load_file('config.yaml')
config$samples <- lapply(config$samples, as.character)
context_list <- c('CpG', 'CHH', 'CHG')

pdf('figures/methylKit_PCA.pdf')

for (context in context_list) {
  cat('check', context, '...\n')
  # read the methylation call files
  file_list <- as.list(paste0('output/', config$samples, '_R1_paired_bismark_bt2_pe_sort.deduplicated.', context, '_report.txt.gz'))
  raw <- methRead(file_list, sample.id=config$samples, assembly=config$assembly, treatment=config$treatment, 
                  context=context, pipeline = 'bismarkCytosineReport', header=F)
  
  # discards bases that have coverage below 10X 
  # and also discards the bases that have more than 
  # 99.9th percentile of coverage in each sample.
  filtered <- filterByCoverage(raw, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
  # merge all samples to one object for base-pair locations that are covered in all samples. 
  meth <- unite(filtered)
  # PCA analysis
  PCASamples(meth)
  
  # write to tsv file
  colnames(meth) <- c(colnames(meth)[1:4], paste0(rep(unlist(config$samples),each=3), '.', c('coverage', 'numCs', 'numTs')))
  meth2 <- getData(meth)[, -grep('coverage', colnames(meth))]
  meth2$start <- meth2$start - 1
  write.table(meth2, paste0('tables/all_', context, '_counts.tsv'), sep='\t', quote = F, row.names = F)
}

dev.off()
