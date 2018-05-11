argv <- commandArgs(T)

load(argv[1])
dmr <- meth[rowSums(meth[,grep('_adj.pvalue', colnames(meth))] < 0.05) > 0, ]
write.table(dmr[,1:3], paste0(argv[2], '.tmp.bed'), quote = F, row.names = F, sep='\t', col.names = F)
system(paste0('bedtools merge -i ', argv[2], '.tmp.bed -d 100 > ', argv[2]))
system(paste0('rm ', argv[2], '.tmp.bed'))
