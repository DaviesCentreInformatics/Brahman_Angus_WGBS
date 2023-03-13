library(methylKit)

treat.file.list <- list('./BxA/F13.Consensus.CpGs.methylKit.gz',
                       './BxA/F62.Consensus.CpGs.methylKit.gz',
                       './BxA/F77.Consensus.CpGs.methylKit.gz',
                       './BxA/F8.Consensus.CpGs.methylKit.gz',
                       './BxA/F80.Consensus.CpGs.methylKit.gz',
                       './BxA/F91.Consensus.CpGs.methylKit.gz')

ctrl.file.list <- list('./AxB/F100.Consensus.CpGs.methylKit.gz',
                        './AxB/F104.Consensus.CpGs.methylKit.gz',
                        './AxB/F106.Consensus.CpGs.methylKit.gz',
                        './AxB/F61.Consensus.CpGs.methylKit.gz',
                        './AxB/F74.Consensus.CpGs.methylKit.gz',
                        './AxB/F97.Consensus.CpGs.methylKit.gz')
################################################################################

file.list <- append(treat.file.list, ctrl.file.list)

############################### Region Analysis ################################
methObj <- methRead(file.list,
                    sample.id = list("F13","F62","F77","F8","F80","F91",
                                     "F100","F104","F106","F61","F74","F97"),
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

meth <- unite(methObj, destrand = TRUE,
             min.per.group = 5L,
             mc.cores = 6)#,
             #save.db = T,
             #suffix = 'min5_destrand')

makeMethylDB(meth, './ABvsBA.consensus.base')
meth <- readMethylDB('./ABvsBA.consensus.base/methylBase_min5_destrand.txt.bgz')

clusterSamples(meth, dist='correlation',method='ward',plot=T)

# Set male = 0, female = 1
covariates <- data.frame(sex=c(1,1,0,1,0,0,
                               1,1,1,0,0,0))

meth_Diff <- calculateDiffMeth(meth, covariates = covariates,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth_bp',
                               slim = T)

meth_tiles = tileMethylCounts(meth,
                              win.size = 1000,
                              step.size = 500,
                              cov.bases = 1,
                              mc.cores = 6,
                              save.db=TRUE,
                              suffix = 'tiled_w1000s500')

# Set male = 0, female = 1
covariates = data.frame(sex=c(1,1,1,0,0,0,
                              1,1,0,1,0,0))

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
write.table(up_dmrs_df, file = './BAvsAB.base/BAvsAB_DMRs_10pHyper.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

down_dmrs_granges <- as(down_dmrs, "GRanges")
down_dmrs_df <- data.frame(seqnames=seqnames(down_dmrs_granges),
                           starts=start(down_dmrs_granges)-1,
                           ends=end(down_dmrs_granges),
                           names=down_dmrs_granges$meth.diff,
                           scores=down_dmrs_granges$qvalue,
                           strands=strand(down_dmrs_granges))
options(scipen = 100, digits=4)
write.table(up_dmrs_df, file = './BAvsAB.base/BAvsAB_DMRs_10pHypo.bed',
            quote = F, sep='\t',row.names = F, col.names = F)
