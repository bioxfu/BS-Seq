library(methylKit)
library(yaml)

config <- yaml.load_file('config.yaml')
context_list <- c('CpG', 'CHG', 'CHH')
win_size <- config$win_size
step_size <- config$step_size

DMR <- list(CpG=list(), CHG=list(), CHH=list())

for (context in context_list) {
  for (i in 1:length(config$versus_0)) {
    control <- config$versus_0[i]
    treat <- config$versus_1[i]
    cat(context, ':', control, 'vs', treat, '....\n')  
    control_sample <- unlist(strsplit(control,',')) 
    treat_sample <- unlist(strsplit(treat,',')) 
    samples <- c(control_sample, treat_sample)
    treatments <- c(rep(0, length(control_sample)), rep(1, length(treat_sample)))
    # read the methylation call files and store them as flat file database
    file_list <- as.list(paste0('output/', samples, '_R1_paired_bismark_bt2_pe_sort.deduplicated.', context, '_report.txt.gz'))
    raw <- methRead(file_list, sample.id=as.list(samples), assembly=config$assembly, treatment=treatments, 
                    context=context, pipeline = 'bismarkCytosineReport', header=F)
    
    # discards bases that have coverage below 10X 
    # and also discards the bases that have more than 
    # 99.9th percentile of coverage in each sample.
    filtered <- filterByCoverage(raw, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
    
    ## DMR
    # summarize methylation information over tiling windows rather than
    # doing base-pair resolution analysis.
    tiles <- tileMethylCounts(filtered, win.size=win_size, step.size=step_size)
    meth_tiles <- unite(tiles)
    
    DMR[[context]][[i]] <- calculateDiffMeth(meth_tiles, mc.cores = 30)
    names(DMR[[context]])[i] <- paste0(control, '_vs_', treat)
    
    ## garbage collection
    rm(list = c('raw', 'filtered', 'tiles', 'meth_tiles'))
    gc()
  }
}

save(list = c('DMR'), file = paste0('RData/methylKit_DMR_w', win_size, '_s', step_size, '.RData'))
   