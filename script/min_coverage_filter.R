library(yaml)

argv <- commandArgs(T)
input <- argv[1]
output <- argv[2]

config <- yaml.load_file('config.yaml')
meth <- read.table(input, header = T)
value <- meth[-(1:2)]
colnames(value) <- paste0(rep(config$groups, each=2), '@', colnames(value))

new_value <- NULL
uniq_group <- unique(config$groups)
for (group in uniq_group) {
  new_value <- cbind(new_value, rowSums(value[, grep(paste0(group,'@'), grep('\\.methyCs', colnames(value), value = T), value = T)]))
  new_value <- cbind(new_value, rowSums(value[, grep(paste0(group,'@'), grep('\\.unmethyCs', colnames(value), value = T), value = T)]))
}
colnames(new_value) <- paste0(rep(uniq_group, each=2), c('.methyCs', '.unmethyCs'))
new_meth <- cbind(meth[1:2], new_value)

coverage <- NULL
for (i in seq(1, ncol(new_value), 2)) {
  coverage <- cbind(coverage, new_value[, i] + new_value[, i+1])
}

filt_idx <- which(rowSums(coverage >= config$min_coverage) == ncol(coverage))
filt_meth <- new_meth[filt_idx,]

write.table(filt_meth, output, sep = '\t', quote = F, row.names = F)
