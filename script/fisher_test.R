treats <- c('70', 'A1', 'A2', 'A3', 'B', 'C')
ctrl <- 'R'

cpg <- read.table('tables/all_CpG_counts.tsv', header = T)
chg <- read.table('tables/all_CHG_counts.tsv', header = T)
chh <- read.table('tables/all_CHH_counts.tsv', header = T)

cpg$context <- 'CpG'
chg$context <- 'CHG'
chh$context <- 'CHH'

meth <- rbind(cpg, chg, chh)

for (treat in treats) {
  cat('testing', treat, '....\n')
  dfm <- data.frame(treat_Cs = rowSums(meth[grep(paste0(treat, '.+numCs'), colnames(meth))]),
                    treat_Ts = rowSums(meth[grep(paste0(treat, '.+numTs'), colnames(meth))]),
                    ctrl_Cs = rowSums(meth[grep(paste0(ctrl, '.+numCs'), colnames(meth))]),
                    ctrl_Ts = rowSums(meth[grep(paste0(ctrl, '.+numTs'), colnames(meth))]))
  mat <- as.matrix(t(dfm))
  
  p <- c()
  for (i in 1:ncol(mat)) {
    p[i] <- fisher.test(matrix(mat[, i],nrow=2))$p.value
  }
  p_adj <- p.adjust(p, method = 'fdr')
  
  meth[paste0(treat,'_vs_',ctrl, '_diff')] <- round(dfm$treat_Cs / (dfm$treat_Cs+dfm$treat_Ts) - dfm$ctrl_Cs / (dfm$ctrl_Cs+dfm$ctrl_Ts),4)
  meth[paste0(treat,'_vs_',ctrl, '_pvalue')] <- round(p, 4)
  meth[paste0(treat,'_vs_',ctrl, '_adj.pvalue')] <- round(p_adj, 4)
}

dmc <- data.frame(pos=paste0(meth$chr, ':', meth$end))
dmc <- cbind(dmc, meth)
dmc$chr <- sub('Chr', '', dmc$chr)
dmc <- dmc[rowSums(dmc[,grep('_pvalue', colnames(dmc))] < 0.05) > 0, ]
colnames(dmc) <- sub('numCs', 'methylatedCs', colnames(dmc))
colnames(dmc) <- sub('numTs', 'unmethylatedCs', colnames(dmc))

write.table(dmc[,1:5], 'tables/peaks.bed', quote = F, row.names = F, sep='\t', col.names = F)

system('annotatePeaks.pl tables/peaks.bed tair10 > tables/peaks.anno')

anno <- read.table('tables/peaks.anno', sep='\t', header = T, row.names = 1, quote = '', comment.char = "")
anno <- anno[, c('Annotation', 'Detailed.Annotation', 'Distance.to.TSS', 'Nearest.PromoterID', 'Gene.Name', 'Gene.Alias', 'Gene.Description', 'Gene.Type')]
total_result_anno <- merge(dmc[-(2:4)], anno, by.x=1, by.y=0)
write.table(total_result_anno, 'tables/DMC_fisher_test_anno.tsv', sep='\t', row.names = F, quote = F)
