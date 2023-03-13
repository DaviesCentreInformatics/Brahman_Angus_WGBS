library(methylKit)
setwd('./aligned2Angus/')

wd <- getwd()
assertthat::are_equal(wd, '/Users/callummacphillamy/PhD/Brah-Ang_methylation/methylation/aligned2Angus')
treat.file.list <- list('./AxB/F100.Consensus.CpGs.methylKit.gz',
                        './AxB/F104.Consensus.CpGs.methylKit.gz',
                        './AxB/F106.Consensus.CpGs.methylKit.gz',
                        './AxB/F61.Consensus.CpGs.methylKit.gz',
                        './AxB/F74.Consensus.CpGs.methylKit.gz',
                        './AxB/F97.Consensus.CpGs.methylKit.gz')

ctrl.file.list <- list('./AxA/F103.Consensus.CpGs.methylKit.gz',
                       './AxA/F105.Consensus.CpGs.methylKit.gz',
                       './AxA/F52.Consensus.CpGs.methylKit.gz',
                       './AxA/F53.Consensus.CpGs.methylKit.gz',
                       './AxA/F60.Consensus.CpGs.methylKit.gz',
                       './AxA/F7.Consensus.CpGs.methylKit.gz')
################################################################################

file.list <- append(treat.file.list, ctrl.file.list)

############################### Region Analysis ################################
methObj <- methRead(file.list,
                    sample.id = list("F100","F104","F106","F61","F74","F97",
                                     "F103","F105","F152","F53","F60","F7"),
                    assembly='UOA-Angus',
                    treatment=c(1,1,1,1,1,1,
                                0,0,0,0,0,0),
                    context='CpG',
                    resolution = 'base')#,
#dbtype="tabix",
#dbdir='./BAvsAB_base')

methObj <- filterByCoverage(methObj,
                            lo.count = 10, lo.perc=NULL,
                            hi.count = NULL, hi.perc = 99.9,
                            save.db = F)

meth <- normalizeCoverage(methObj)

meth = unite(meth, destrand = TRUE,
             min.per.group = 5L,
             mc.cores = 6)#,
#save.db = T,
#suffix = 'min5_destrand')

makeMethylDB(meth, './AAvsAB.consensus.base')
meth <- readMethylDB('./AAvsAB.consensus.base/methylBase_min5_destrand.txt.bgz')

#clusterSamples(meth, dist='correlation',method='ward',plot=T)

# Set male = 0, female = 1
covariates <- data.frame(sex=c(1,1,1,0,0,0,
                               0,0,0,1,1,1))

meth_Diff <- calculateDiffMeth(meth, covariates = covariates,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth_bp',
                               slim = T)

meth_tiles = tileMethylCounts(meth,
                              mc.cores = 6,
                              save.db=TRUE,
                              suffix = 'tiled')

# Set male = 0, female = 1

covariates = data.frame(sex=c(1,1,1,0,0,0,
                              0,0,0,1,1,1))

meth_Diff <- calculateDiffMeth(meth_tiles, covariates = covariates,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth')



up_dmrs <- getMethylDiff(meth_Diff, difference = 10, type='hyper',save.db = F)
down_dmrs <- getMethylDiff(meth_Diff, difference = 10, type='hypo',save.db = F)

up_dmrs_granges <- as(up_dmrs, "GRanges")
up_dmrs_df <- data.frame(seqnames=seqnames(up_dmrs_granges),
                         starts=start(up_dmrs_granges)-1,
                         ends=end(up_dmrs_granges),
                         names=up_dmrs_granges$meth.diff,
                         scores=up_dmrs_granges$qvalue,
                         strands=strand(up_dmrs_granges))
options(scipen = 100, digits=4)
write.table(up_dmrs_df, file = './ABvsAA.base/ABvsAA_DMRs_10pHyper.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

down_dmrs_granges <- as(down_dmrs, "GRanges")
down_dmrs_df <- data.frame(seqnames=seqnames(down_dmrs_granges),
                           starts=start(down_dmrs_granges)-1,
                           ends=end(down_dmrs_granges),
                           names=down_dmrs_granges$meth.diff,
                           scores=down_dmrs_granges$qvalue,
                           strands=strand(down_dmrs_granges))
options(scipen = 100, digits=4)
write.table(up_dmrs_df, file = './ABvsAA.base/ABvsAA_DMRs_10pHypo.bed',
            quote = F, sep='\t',row.names = F, col.names = F)
