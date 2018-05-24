library(yaml)
#library(BiocParallel)

config <- yaml.load_file('config.yaml')
treats <- strsplit(config$fisher_test$treat, ',')[[1]]
ctrls <- strsplit(config$fisher_test$ctrl, ',')[[1]]

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

meth <- read.table(input, header = T)
meth2 <- meth[1:2]

meth_test <- list()
for (i in 1:length(treats)) {
  treat <- treats[i]
  ctrl <- ctrls[i]
  cat('testing', treat, ' vs ', ctrl, '....\n')
  dfm <- data.frame(treat_Cs = meth[grep(paste0(treat, '\\.methyCs'), colnames(meth), value=T)],
                    treat_Ts = meth[grep(paste0(treat, '\\.unmethyCs'), colnames(meth), value=T)],
                    ctrl_Cs = meth[grep(paste0(ctrl, '\\.methyCs'), colnames(meth), value=T)],
                    ctrl_Ts = meth[grep(paste0(ctrl, '\\.unmethyCs'), colnames(meth), value=T)])
  colnames(dfm) <- c('treat_Cs', 'treat_Ts', 'ctrl_Cs', 'ctrl_Ts')
  treat_level <- dfm$treat_Cs / (dfm$treat_Cs + dfm$treat_Ts) + 0.0001
  ctrl_level <- dfm$ctrl_Cs / (dfm$ctrl_Cs + dfm$ctrl_Ts) + 0.0001
  FC <- treat_level / ctrl_level
  log2FC <- log2(FC)

  #p <- bplapply(dfm, function(x){fisher.test(matrix(x,nrow=2))$p.value}, BPPARAM = MulticoreParam(workers = 20, progressbar = TRUE))
  p <- lapply(as.data.frame(t(dfm)), function(x){fisher.test(matrix(x,nrow=2))$p.value})
  p <- unlist(p)
  #p_adj <- p.adjust(p, method = 'fdr')
  
  meth3 <- meth2
  meth3[paste0(treat,'_vs_',ctrl, '.FC')] <- round(FC, 4)
  meth3[paste0(treat,'_vs_',ctrl, '.log2FC')] <- round(log2FC, 4)
  meth3[paste0(treat,'_vs_',ctrl, '.p')] <- round(p, 4)
  #meth3[paste0(treat,'_vs_',ctrl, '.FDR')] <- round(p_adj, 4)
  meth_test[[i]] <- meth3
  names(meth_test)[i] <- paste0(treat, '_vs_', ctrl)
}

save(meth_test, file = output)
