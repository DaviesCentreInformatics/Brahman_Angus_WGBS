library('methylKit')
library('GenomicRanges')


treat.file.list <- list('./AxB/methylKit/F100.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F104.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F106.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F61.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F74.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F97.Consensus.CpGs.methylKit.gz')

ctrl.file.list <- list('./BxB/methylKit/F22.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F46.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F56.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F65.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F78.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F99.Consensus.CpGs.methylKit.gz')

file.list <- append(treat.file.list, ctrl.file.list)

############################### BP-DMR Analysis ################################
methObj <- methRead(file.list,
                    sample.id = list("F100","F104","F106","F61","F74","F97",
                                     "F22","F46","F56","F65","F78","F99"),
                    assembly='UOA-Brahman',
                    treatment=c(1,1,1,1,1,1,
                                0,0,0,0,0,0),
                    context='CpG',
                    resolution = 'base')#,
                    #dbtype="tabix",
                    #dbdir='./BBvsAB.consensus_base')

meth <- filterByCoverage(methObj,
                         lo.count = 10, lo.perc=NULL,
                         hi.count = NULL, hi.perc = 99.9,
                         save.db = F)

meth <- normalizeCoverage(meth)

meth <- unite(meth, destrand = T,
             min.per.group = 5L,
             mc.cores = 6)#,
             #save.db = T,
             #suffix = 'min5_destrand')

makeMethylDB(meth, dbdir='./BBvsAB.consensus.base/')

meth <- readMethylDB('./BBvsAB.base/methylBase_min5_destrand.txt.bgz')

covariates <- data.frame(sex=c(1,1,1,0,0,0,
                               1,0,0,1,0,1))

meth_Diff <- calculateDiffMeth(meth, covariates = covariates,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth_bp',
                               slim = T)

meth_tiles <- tileMethylCounts(meth,
                              mc.cores = 6,
                              save.db=TRUE,
                              suffix = 'tiled')

# Set male = 0, female = 1
covariates <- data.frame(sex=c(1,1,1,0,0,0,
                              1,0,0,1,0,1))

meth_Diff <- calculateDiffMeth(meth_tiles, covariates = covariates,
                              mc.cores = 6,
                              save.db = T,
                              suffix = 'diffMeth',
                              slim = T)

meth_Diff <- readMethylDB('./BBvsAB.base/methylDiff_min5_destrand_tiled_w1000s500_diffMeth.txt.bgz')

up_dmrs <- getMethylDiff(meth_Diff, difference=10, type='hyper', save.db = F)
down_dmrs <- getMethylDiff(meth_Diff, difference = 10, type='hypo', save.db=F)

up_dmrs_granges <- as(up_dmrs, "GRanges")
up_dmrs_df <- data.frame(seqnames=seqnames(up_dmrs_granges),
                         starts=start(up_dmrs_granges)-1,
                         ends=end(up_dmrs_granges),
                         names=up_dmrs_granges$meth.diff,
                         scores=up_dmrs_granges$qvalue,
                         strands=strand(up_dmrs_granges))
options(scipen = 100, digits=4)
write.table(up_dmrs_df, file = './BBvsAB.base/BBvsAB_DMRs_10pHyper.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

down_dmrs_granges <- as(down_dmrs, "GRanges")
down_dmrs_df <- data.frame(seqnames=seqnames(down_dmrs_granges),
                           starts=start(down_dmrs_granges)-1,
                           ends=end(down_dmrs_granges),
                           names=down_dmrs_granges$meth.diff,
                           scores=down_dmrs_granges$qvalue,
                           strands=strand(down_dmrs_granges))
options(scipen = 100, digits=4)
write.table(down_dmrs_df, file = './BBvsAB.base/BBvsAB_DMRs_10pHypo.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

