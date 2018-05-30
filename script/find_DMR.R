library(yaml)
config <- yaml.load_file('config.yaml')

cluster <- function(m) {
  x <- m$pos
  steps <- diff(x)
  b <- which(steps > 100)
  
  start_idx <- c(1, b + 1)
  end_idx <- c(b, length(x))
  start_pos <- x[start_idx]
  end_pos <- x[end_idx]
  
  id <- paste0(unique(m$chrom), ':', start_pos, '-', end_pos)
  m$window <- rep(id, end_idx - start_idx + 1)
  return(m)  
}

contexts <- c('CpG', 'CHH', 'CHG')

for (context in contexts) {
  
  load(paste0('fisher_test/methyC_', context, '.RData'))
  
  for (x in 1:length(meth_test_all)) {
    cat(paste(context, names(meth_test_all)[x], '....\n'))
    meth <- meth_test_all[[x]]
    meth$FDR <- p.adjust(meth[, 5])
    meth_chrom <- split(meth, meth$chrom)
    
    meth2 <- NULL
    for (i in 1:length(meth_chrom)) {
      m <- meth_chrom[[i]]
      m <- m[order(m$pos), ]
      m_cluster <- cluster(m)
      m_cluster$cutoff <- m_cluster[, 6] < 0.05 & abs(m_cluster[, 4]) > log2(config$fold_change)
      filt <- tapply(m_cluster$cutoff, m_cluster$window, sum)
    
      DMC <- m_cluster[m_cluster$window %in% names(filt[filt >= 5]), ]
      if (nrow(DMC) > 0) {
        DMR <- aggregate(DMC[, 3], list(window=DMC$window), mean)
        y <- merge(DMC[DMC$cutoff==TRUE,], DMR)
        DMR_start <- tapply(y$pos, y$window, function(x){min(x)})
        DMR_end <- tapply(y$pos, y$window, function(x){max(x)})
        DMC_num <- tapply(y$pos, y$window, function(x){length(x)})
        new_bound <- cbind(names(meth_chrom)[i], DMR_start, DMR_end, DMC_num)
        meth2 <- rbind(meth2, merge(new_bound, DMR, by.x=0, by.y=1))
      }
    }
    if (!is.null(meth2)) {
      meth2$log2FC <- round(log2(meth2$x), 4)
      meth2 <- meth2[abs(meth2$log2FC) > log2(config$fold_change), -c(1,6)]
      colnames(meth2)[1] <- 'chrom'
      meth2$DMR_size <- as.numeric(as.vector(meth2$DMR_end)) - as.numeric(as.vector(meth2$DMR_start)) + 1
      meth2$direction <- meth2$log2FC
      meth2$direction[meth2$log2FC > 0] <- 'up'
      meth2$direction[meth2$log2FC < 0] <- 'down'
    }
    write.table(meth2, paste0('tables/', context, '_DMR_', gsub('_+', '.', names(meth_test_all)[x]), '.tsv'), row.names = F, sep='\t', quote = F)
  }
}

