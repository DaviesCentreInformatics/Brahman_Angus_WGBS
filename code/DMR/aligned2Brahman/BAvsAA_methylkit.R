library(methylKit)
treat.file.list <- list('./AxA/methylKit/F103_1.sorted.markDups_CpG.methylKit',
                        './AxA/methylKit/F105_1.sorted.markDups_CpG.methylKit',
                        './AxA/methylKit/F52_1.sorted.markDups_CpG.methylKit',
                        './AxA/methylKit/F53_1.sorted.markDups_CpG.methylKit',
                        './AxA/methylKit/F60_1.sorted.markDups_CpG.methylKit',
                        './AxA/methylKit/F7_1.sorted.markDups_CpG.methylKit')

ctrl.file.list <- list('./BxA/methylKit/F13_1.sorted.markDups_CpG.methylKit',
                       './BxA/methylKit/F62_1.sorted.markDups_CpG.methylKit',
                       './BxA/methylKit/F77_1.sorted.markDups_CpG.methylKit',
                       './BxA/methylKit/F8_1.sorted.markDups_CpG.methylKit',
                       './BxA/methylKit/F80_1.sorted.markDups_CpG.methylKit',
                       './BxA/methylKit/F91_1.sorted.markDups_CpG.methylKit')
################################################################################

file.list <- append(treat.file.list, ctrl.file.list)

############################### Region Analysis ################################
methObj <- methRead(file.list,
                    sample.id = list("F103","F105","F152","F53","F60","F7",
                                     "F13","F62","F77","F8","F80","F91"),
                    assembly='UOA-Brahman',
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

makeMethylDB(meth, './BAvsAA.base')
meth <- readMethylDB('./BAvsAA.base/methylBase_min5_destrand_tiled.txt.bgz')

#clusterSamples(meth, dist='correlation',method='ward',plot=T)

# Set male = 0, female = 1
covariates <- data.frame(sex=c(0,0,0,1,1,1,
                               1,1,0,1,0,0))

meth_Diff <- calculateDiffMeth(meth,
                               covariates = covariates,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth_bp',
                               slim = T)

meth_tiles = tileMethylCounts(meth,
                              mc.cores = 6,
                              save.db=TRUE,
                              suffix = 'tiled')

# Set male = 0, female = 1

covariates = data.frame(sex=c(0,0,0,1,1,1,
                              1,1,0,1,0,0))

meth_Diff <- calculateDiffMeth(meth, covariates = covariates,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth')

meth_Diff <- readMethylDB('./BAvsAA.base/methylDiff_min5_destrand_tiled_w1000s500_diffMeth.txt.bgz')

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
write.table(up_dmrs_df, file = './BAvsAA.base/BAvsAA_DMRs_10pHyper.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

down_dmrs_granges <- as(down_dmrs, "GRanges")
down_dmrs_df <- data.frame(seqnames=seqnames(down_dmrs_granges),
                           starts=start(down_dmrs_granges)-1,
                           ends=end(down_dmrs_granges),
                           names=down_dmrs_granges$meth.diff,
                           scores=down_dmrs_granges$qvalue,
                           strands=strand(down_dmrs_granges))
options(scipen = 100, digits=4)
write.table(up_dmrs_df, file = './BAvsAA.base/BAvsAA_DMRs_10pHypo.bed',
            quote = F, sep='\t',row.names = F, col.names = F)
