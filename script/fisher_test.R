library(yaml)

config <- yaml.load_file('config.yaml')
treats <- strsplit(config$fisher_test$treat, ',')[[1]]
ctrls <- strsplit(config$fisher_test$ctrl, ',')[[1]]

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

meth <- read.table(input, header = T)

for (i in 1:length(treats)) {
  treat <- treats[i]
  ctrl <- ctrls[i]
  cat('testing', treat, ' vs ', ctrl, '....\n')
  dfm <- data.frame(treat_Cs = rowSums(meth[grep(paste0(treat, '.+methyCs'), colnames(meth))]),
                    treat_Ts = rowSums(meth[grep(paste0(treat, '.+unmethyCs'), colnames(meth))]),
                    ctrl_Cs = rowSums(meth[grep(paste0(ctrl, '.+methyCs'), colnames(meth))]),
                    ctrl_Ts = rowSums(meth[grep(paste0(ctrl, '.+unmethyCs'), colnames(meth))]))
  mat <- as.matrix(t(dfm))
  
  p <- c()
  total <- ncol(mat)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for (i in 1:total) {
    setTxtProgressBar(pb, i)
    p[i] <- fisher.test(matrix(mat[, i],nrow=2))$p.value
  }
  p_adj <- p.adjust(p, method = 'fdr')
  
  meth[paste0(treat,'_vs_',ctrl, '_diff')] <- round(dfm$treat_Cs / (dfm$treat_Cs+dfm$treat_Ts) - dfm$ctrl_Cs / (dfm$ctrl_Cs+dfm$ctrl_Ts),4)
  meth[paste0(treat,'_vs_',ctrl, '_pvalue')] <- round(p, 4)
  meth[paste0(treat,'_vs_',ctrl, '_adj.pvalue')] <- round(p_adj, 4)
}

save(meth, file = output)
