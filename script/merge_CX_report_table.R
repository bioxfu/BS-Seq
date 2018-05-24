library(yaml)

config <- yaml.load_file('config.yaml')
config$samples <- lapply(config$samples, as.character)
context_list <- c('CpG', 'CHH', 'CHG')

for (context in context_list) {
  files <- paste0('output/', config$samples, '_bismark_sort.deduplicated.', context, '_report.txt.gz')
  count <- do.call("cbind", lapply(files, read.table, sep='\t'))
  count <- count[c(1:2, grep('V[45]', colnames(count)))]  
  colnames(count) <- c('chrom', 'pos', paste0(rep(config$samples, each=2), c('.methyCs', '.unmethyCs')))
  count$chrom <- sub('Chr', '', count$chrom)
  write.table(count, paste0('tables/methyC_', context, '_counts.tsv'), sep='\t', quote = F, row.names = F)
}
