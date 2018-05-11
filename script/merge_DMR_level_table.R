library(yaml)

config <- yaml.load_file('config.yaml')
config$samples <- lapply(config$samples, as.character)
context_list <- c('CpG', 'CHH', 'CHG')

for (context in context_list) {
  files <- paste0('tables/DMR/methyC_windows_200_step_50_', context, '_DMR_', config$samples, '.bed')
  count <- do.call("cbind", lapply(files, read.table, sep='\t'))
  count <- count[sort(c(1,4,5,6, seq(7,ncol(count),8), seq(8,ncol(count),8)))]  
  colnames(count) <- c('chrom', 'start', 'end', 'num_C', paste0(rep(config$samples, each=2), c('.methyCs', '.unmethyCs')))
  write.table(count, paste0('tables/methyC_windows_200_step_50_', context, '_DMR_counts.tsv'), sep='\t', quote = F, row.names = F)
}

