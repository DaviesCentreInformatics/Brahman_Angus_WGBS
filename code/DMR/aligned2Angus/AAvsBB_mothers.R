library('methylKit')

getwd()
setwd('/Users/callummacphillamy/PhD/Brah-Ang_methylation/methylation/aligned2Angus/')
getwd()
# Angus mothers
ctrl.file.list <- list('./AxA/F103.Consensus.CpGs.methylKit.gz',
                        './AxA/F105.Consensus.CpGs.methylKit.gz',
                        './AxA/F52.Consensus.CpGs.methylKit.gz',
                        './AxA/F53.Consensus.CpGs.methylKit.gz',
                        './AxA/F60.Consensus.CpGs.methylKit.gz',
                        './AxA/F7.Consensus.CpGs.methylKit.gz',
                        './BxA/F13.Consensus.CpGs.methylKit.gz',
                        './BxA/F62.Consensus.CpGs.methylKit.gz',
                        './BxA/F77.Consensus.CpGs.methylKit.gz',
                        './BxA/F8.Consensus.CpGs.methylKit.gz',
                        './BxA/F80.Consensus.CpGs.methylKit.gz',
                        './BxA/F91.Consensus.CpGs.methylKit.gz')

# Brahman mothers
treat.file.list <- list('./BxB/F22.Consensus.CpGs.methylKit.gz',
                       './BxB/F46.Consensus.CpGs.methylKit.gz',
                       './BxB/F56.Consensus.CpGs.methylKit.gz',
                       './BxB/F65.Consensus.CpGs.methylKit.gz',
                       './BxB/F78.Consensus.CpGs.methylKit.gz',
                       './BxB/F99.Consensus.CpGs.methylKit.gz',
                       './AxB/F100.Consensus.CpGs.methylKit.gz',
                       './AxB/F104.Consensus.CpGs.methylKit.gz',
                       './AxB/F106.Consensus.CpGs.methylKit.gz',
                       './AxB/F61.Consensus.CpGs.methylKit.gz',
                       './AxB/F74.Consensus.CpGs.methylKit.gz',
                       './AxB/F97.Consensus.CpGs.methylKit.gz')

file.list <- append(treat.file.list, ctrl.file.list)

############################### BP-DMR Analysis ################################
methObj <- methRead(file.list,
                    sample.id = list("F22","F46","F56","F65","F78","F99","F100","F104","F106","F61","F74","F97",
                                     "F103","F105","F52","F53","F60","F7",'F13','F62',"F77","F8","F80","F91"),
                    assembly='UOA-Angus',
                    treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,
                                0,0,0,0,0,0,0,0,0,0,0,0),
                    context='CpG',
                    resolution = 'base')#,
                    #dbtype="tabix",
                    #dbdir='./AAvsBB_mothers_base')

meth <- filterByCoverage(methObj,
                         lo.count = 10, lo.perc=NULL,
                         hi.count = NULL, hi.perc = 99.9)#,
#save.db = TRUE)

meth <- normalizeCoverage(meth)

meth <- unite(meth, destrand = T,
              min.per.group = 10L,
              mc.cores = 6)#,
#save.db = T,
#suffix = 'min5_destrand')

makeMethylDB(meth,'./AAvsBB_mothers_base.consensus')

meth <- readMethylDB('./AAvsBB_mothers_base.consensus/methylBase_min10_destrand.txt.bgz')

# Set male = 0, female = 1
covariates <- data.frame(sex=c(1,0,0,1,0,1,1,1,1,0,0,0,
                               0,0,0,1,1,1,1,1,0,1,0,0))

meth_Diff <- calculateDiffMeth(meth,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth_bp',
                               slim = T)

meth_tiles <- tileMethylCounts(meth,
                               mc.cores = 6,
                               save.db=TRUE,
                               suffix = 'tiled')


meth_Diff <- calculateDiffMeth(meth_tiles,
                               covariates=covariates,
                               mc.cores = 10,
                               save.db = T,
                               suffix = 'diffMeth',
                               slim = T)

meth <- readMethylDB('./aligned2Angus/AAvsBB_mothers_base.consensus/methylDiff_min10_destrand_tiled_diffMeth.txt.bgz')
dmrs <- getMethylDiff(meth_Diff, difference = 10, save.db = F, type='hypo')
dim(dmrs)
dmrs
