library(methylKit)
library(dplyr)

load('RData/methylKit_DMR.RData')

context_list <- c('CpG', 'CHG', 'CHH')
total_result <- NULL

for (context in context_list) {
  for (i in 1:length(DMR[[context]])) {

    DMR.hyper <- getMethylDiff(DMR[[context]][[i]], difference=10, qvalue=0.05, type="hyper")
    DMR.hypo <- getMethylDiff(DMR[[context]][[i]], difference=10, qvalue=0.05, type="hypo")

    result1 <- data.frame(window = paste0(DMR.hyper$chr, ':', DMR.hyper$start, '-', DMR.hyper$end),
                          chrom = sub('Chr', '', DMR.hyper$chr),
                          start = DMR.hyper$start,
                          end = DMR.hyper$end,
                          type = context,
                          change = 'hyper',
                          VS = names(DMR[[context]])[i])
    
    result2 <- data.frame(window = paste0(DMR.hypo$chr, ':', DMR.hypo$start, '-', DMR.hypo$end),
                          chrom = sub('Chr', '', DMR.hypo$chr),
                          start = DMR.hypo$start,
                          end = DMR.hypo$end,
                          type = context,
                          change = 'hypo',
                          VS = names(DMR[[context]])[i])
    
    total_result <- rbind(total_result, rbind(result1, result2))

  }
}

total_result_summary <- group_by(total_result, window, chrom, start, end, type, change) %>% summarise(vs = paste(VS, collapse = '|'))

peak = unique(total_result_summary[, 1:4])
peak$strand <- '+'

write.table(peak, 'tables/peaks.bed', sep='\t', row.names = F, col.names = F, quote = F)

system('annotatePeaks.pl tables/peaks.bed tair10 > tables/peaks.anno')

anno <- read.table('tables/peaks.anno', sep='\t', header = T, row.names = 1, quote = '', comment.char = "")
anno <- anno[, c('Annotation', 'Detailed.Annotation', 'Distance.to.TSS', 'Nearest.PromoterID', 'Gene.Name', 'Gene.Alias', 'Gene.Description', 'Gene.Type')]
total_result_anno <- merge(total_result_summary, anno, by.x=1, by.y=0)

write.table(total_result_anno, 'tables/methylKit_DMR_anno.tsv', sep='\t', row.names = F, quote = F)
