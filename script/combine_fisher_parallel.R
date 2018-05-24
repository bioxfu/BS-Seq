library(yaml)
config <- yaml.load_file('config.yaml')
contexts <- c('CpG', 'CHH', 'CHG')

for (context in contexts) {
  Rdata <- dir(path = paste0('fisher_test/methyC_', context, '/'), pattern = '.RData', full.names = T)
  load(Rdata[1])
  meth_test_all <- meth_test
  
  for (x in 2:length(Rdata)) {
    cat(Rdata[x], '...\n')
    load(Rdata[x])
    for (i in 1:length(meth_test)) {
      meth_test_all[[i]] <- rbind(meth_test_all[[i]], meth_test[[i]])
    }
  }
  save(meth_test_all, file = paste0('fisher_test/methyC_', context, '.RData'))
}

