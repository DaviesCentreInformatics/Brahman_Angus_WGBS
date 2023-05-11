library('methylKit')
setwd('/Users/callummacphillamy/PhD/Brah-Ang_methylation/methylation')
# Angus fathers
treat.file.list <- list('./AxA/methylKit/F103.Consensus.CpGs.methylKit.gz',
                        './AxA/methylKit/F105.Consensus.CpGs.methylKit.gz',
                        './AxA/methylKit/F52.Consensus.CpGs.methylKit.gz',
                        './AxA/methylKit/F53.Consensus.CpGs.methylKit.gz',
                        './AxA/methylKit/F60.Consensus.CpGs.methylKit.gz',
                        './AxA/methylKit/F7.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F100.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F104.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F106.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F61.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F74.Consensus.CpGs.methylKit.gz',
                        './AxB/methylKit/F97.Consensus.CpGs.methylKit.gz')

# Brahman fathers
ctrl.file.list <- list('./BxB/methylKit/F22.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F46.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F56.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F65.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F78.Consensus.CpGs.methylKit.gz',
                       './BxB/methylKit/F99.Consensus.CpGs.methylKit.gz',
                       './BxA/methylKit/F13.Consensus.CpGs.methylKit.gz',
                       './BxA/methylKit/F62.Consensus.CpGs.methylKit.gz',
                       './BxA/methylKit/F77.Consensus.CpGs.methylKit.gz',
                       './BxA/methylKit/F8.Consensus.CpGs.methylKit.gz',
                       './BxA/methylKit/F80.Consensus.CpGs.methylKit.gz',
                       './BxA/methylKit/F91.Consensus.CpGs.methylKit.gz')

file.list <- append(treat.file.list, ctrl.file.list)

############################### BP-DMR Analysis ################################
methObj <- methRead(file.list,
                    sample.id = list("F103","F105","F52","F53","F60","F7","F100","F104","F106","F61","F74","F97",
                                     "F22","F46","F56","F65","F78","F99",'F13','F62',"F77","F8","F80","F91"),
                    assembly='UOA-Brahman',
                    treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,
                                0,0,0,0,0,0,0,0,0,0,0,0),
                    context='CpG',
                    resolution = 'base')#,
                    #dbtype="tabix",
                    #dbdir='./BBvsAA_fathers_base.consensus')

meth <- filterByCoverage(methObj,
                         lo.count = 10, lo.perc=NULL,
                         hi.count = NULL, hi.perc = 99.9)#,
#save.db = TRUE)

meth <- normalizeCoverage(meth)

meth <- unite(meth, destrand = T,
              min.per.group = 10L,
              mc.cores = 10)#,
              #save.db = T,
              #suffix = 'min10_destrand')

makeMethylDB(meth,'./BBvsAA_fathers_base.consensus')

meth <- readMethylDB('./BBvsAA_fathers_base.consensus/methylBase_min10_destrand_tiled.txt.bgz')

# Set male = 0, female = 1
covariates <- data.frame(sex=c(0,0,0,1,1,1,1,1,1,0,0,0,
                               1,0,0,1,0,1,1,1,0,1,0,0))

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
                               covariates = covariates,
                               mc.cores = 10,
                               save.db = T,
                               suffix = 'diffMeth',
                               slim = T)
