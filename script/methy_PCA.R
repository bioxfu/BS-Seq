library(RColorBrewer)
library(yaml)

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]
#input <- 'tables/methyC_CpG_counts.tsv'
#output <- 'figures/methyC_CpG_counts_PCA.pdf'

config <- yaml.load_file('config.yaml')
meth <- read.table(input, header = T)
value <- meth[-(1:2)]

idx <- seq(1, ncol(value), 2)
meth_cover <- matrix(nrow = nrow(value), ncol = length(idx))
for (i in 1:length(idx)) {
  meth_cover[,i] <- value[, idx[i]] + value[, idx[i]+1]
}
rownames(meth_cover) <- rownames(value)
colnames(meth_cover) <- config$samples
filt <- rowSums(meth_cover > config$min_coverage) == ncol(meth_cover)

value <- value[filt, ] + 1
meth_ratio <- matrix(nrow = nrow(value), ncol = length(idx))
for (i in 1:length(idx)) {
  meth_ratio[,i] <- value[, idx[i]] / (value[, idx[i]] + value[, idx[i]+1])
}
rownames(meth_ratio) <- rownames(value)
colnames(meth_ratio) <- config$samples

pca <- prcomp(t(meth_ratio), center = T, scale. = T)

x <- pca$x
prop_var <- round(summary(pca)$importance[2,1:2]*100,0)

set2_cols <- brewer.pal(12, 'Paired')
cols <- rep(set2_cols, each=config$seq_info$replicate)

pdf(output, hei=7, wid=7)
layout(matrix(c(1,2),nrow=1), wid=c(5, 5))
par(mar=c(5,4,4,0))
par(xpd=TRUE)
plot(x[,1], x[,2], xlim=range(x[,1])*1.1, ylim=range(x[,2])*1.1, col=cols, cex=1, type='n',
     xlab=paste0('PC1 (',prop_var[1],'% of Variance)'),
     ylab=paste0('PC2 (',prop_var[2],'% of Variance)')
)
text(x[,1], x[,2], colnames(meth_ratio), col=cols)
par(mar=c(5,0,4,0))
plot.new()
legend('left', unique(config$groups), pch = 1, col=set2_cols, bty='n')
dev.off()

