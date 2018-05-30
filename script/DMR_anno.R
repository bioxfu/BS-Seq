argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

meth <- read.table(input, header = T)
dmr <- data.frame(pos=paste0(meth$chrom, ':', meth$DMR_start, '-', meth$DMR_end))
dmr <- cbind(dmr, meth)

bed <- dmr[1:4]
bed$strand <- '+'
write.table(bed, paste0(input, '.tmp.peaks.bed'), quote = F, row.names = F, sep='\t', col.names = F)
system(paste0('annotatePeaks.pl ', input, '.tmp.peaks.bed tair10 > ', input, '.tmp.peaks.anno'))

anno <- read.table(paste0(input, '.tmp.peaks.anno'), sep='\t', header = T, row.names = 1, quote = '', comment.char = "")
anno <- anno[, c('Annotation', 'Detailed.Annotation', 'Distance.to.TSS', 'Nearest.PromoterID', 'Gene.Name', 'Gene.Alias', 'Gene.Description', 'Gene.Type')]
total_result_anno <- merge(dmr[-(2:4)], anno, by.x=1, by.y=0)
write.table(total_result_anno, output, sep='\t', row.names = F, quote = F)

system(paste0('rm ', input, '.tmp.peaks.*'))

