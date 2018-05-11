library(RColorBrewer)
cols <- brewer.pal(3, 'Set1')

argv <- commandArgs(T)
input <- argv[1]
output1 <- argv[2]
output2 <- argv[3]

load(input)

meth <- meth[rowSums(meth[,grep('_adj.pvalue', colnames(meth))] < 0.05) > 0, ]
dmr <- data.frame(pos=paste0(meth$chrom, ':', meth$start+1, '-', meth$end))
dmr <- cbind(dmr, meth)
dmr$chrom <- sub('Chr', '', dmr$chrom)

bed <- dmr[1:4]
bed$strand <- '+'
write.table(bed, paste0(input, '.tmp.peaks.bed'), quote = F, row.names = F, sep='\t', col.names = F)
system(paste0('annotatePeaks.pl ', input, '.tmp.peaks.bed tair10 > ', input, '.tmp.peaks.anno'))

anno <- read.table(paste0(input, '.tmp.peaks.anno'), sep='\t', header = T, row.names = 1, quote = '', comment.char = "")
anno <- anno[, c('Annotation', 'Detailed.Annotation', 'Distance.to.TSS', 'Nearest.PromoterID', 'Gene.Name', 'Gene.Alias', 'Gene.Description', 'Gene.Type')]
total_result_anno <- merge(dmr[-(2:4)], anno, by.x=1, by.y=0)
write.table(total_result_anno, output1, sep='\t', row.names = F, quote = F)

stat <- data.frame(up = colSums(total_result_anno[,grep('_adj.pvalue', colnames(total_result_anno))] < 0.05 & total_result_anno[,grep('_diff', colnames(total_result_anno))] > 0),
                   down = colSums(total_result_anno[,grep('_adj.pvalue', colnames(total_result_anno))] < 0.05 & total_result_anno[,grep('_diff', colnames(total_result_anno))] < 0))

rownames(stat) <- sub('_adj.pvalue', '', rownames(stat))
stat <- t(as.matrix(stat))

pdf(output2)
bp <- barplot(stat, beside = T, ylab = 'Number of DMR', col = cols[1:2], border = 'white', las=3, ylim = c(0, max(stat)* 1.1))
legend('topleft', rownames(stat), fill = cols[1:2], bty='n', border = F)
text(bp, stat, stat, pos = 3)
dev.off()

system(paste0('rm ', input, '.tmp.peaks.*'))

