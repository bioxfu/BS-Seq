library(methylKit)
library(yaml)

config <- yaml.load_file('config.yaml')
context_list <- c('CpG', 'CHG', 'CHH')

DMC <- list(CpG=list(), CHG=list(), CHH=list())
DMR <- list(CpG=list(), CHG=list(), CHH=list())

for (context in context_list) {
  for (i in 1:length(config$versus_0)) {
    control <- config$versus_0[i]
    treat <- config$versus_1[i]
    samples <- c(control, treat)
    cat(context, ':', control, 'vs', treat, '....\n')  
  
    # read the methylation call files and store them as flat file database
    file_list <- as.list(paste0('output/', samples, '_clean_bismark_bt2.deduplicated.', context, '_report.txt.gz'))
    raw <- methRead(file_list, sample.id=as.list(samples), assembly=config$assembly, treatment=c(0, 1), 
                    context=context, pipeline = 'bismarkCytosineReport', header=F)
    
    # discards bases that have coverage below 10X 
    # and also discards the bases that have more than 
    # 99.9th percentile of coverage in each sample.
    filtered <- filterByCoverage(raw, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
    
    # merge all samples to one object for base-pair locations that are covered in all samples. 
    meth <- unite(filtered)
    
    ## DMC
    DMC[[context]][[i]] <- calculateDiffMeth(meth, num.cores = 20)
    # DMC.hyper <- getMethylDiff(DMC, difference=25, qvalue=0.05, type="hyper")
    # DMC.hypo <- getMethylDiff(DMC, difference=25, qvalue=0.05, type="hypo")
    # DMC.all <- getMethylDiff(DMC, difference=25, qvalue=0.05)
    
    ## DMR
    # summarize methylation information over tiling windows rather than
    # doing base-pair resolution analysis.
    tiles <- tileMethylCounts(filtered, win.size=1000, step.size=1000)
    meth_tiles <- unite(tiles)
    
    DMR[[context]][[i]] <- calculateDiffMeth(meth_tiles, num.cores = 20)
    # DMR.hyper <- getMethylDiff(DMR, difference=25, qvalue=0.05, type="hyper")
    # DMR.hypo <- getMethylDiff(DMR, difference=25, qvalue=0.05, type="hypo")
    # DMR.all <- getMethylDiff(DMR, difference=25, qvalue=0.05)
    
    names(DMC[[context]])[i] <- paste0(control, '_vs_', treat)
    names(DMR[[context]])[i] <- paste0(control, '_vs_', treat)
  }
  
}

save(list = c('DMC', 'DMR'), file = 'RData/methylKit_DMR.RData')

