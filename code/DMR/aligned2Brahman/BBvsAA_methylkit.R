library('methylKit')
#treat.file.list <- list('./AxA/methylKit/F103.Consensus.CpGs.methylKit.gz',
#                        './AxA/methylKit/F105.Consensus.CpGs.methylKit.gz',
#                        './AxA/methylKit/F52.Consensus.CpGs.methylKit.gz',
#                        './AxA/methylKit/F53.Consensus.CpGs.methylKit.gz',
#                        './AxA/methylKit/F60.Consensus.CpGs.methylKit.gz',
#                        './AxA/methylKit/F7.Consensus.CpGs.methylKit.gz')

#ctrl.file.list <- list('./BxB/methylKit/F22.Consensus.CpGs.methylKit.gz',
#                       './BxB/methylKit/F46.Consensus.CpGs.methylKit.gz',
#                       './BxB/methylKit/F56.Consensus.CpGs.methylKit.gz',
#                       './BxB/methylKit/F65.Consensus.CpGs.methylKit.gz',
#                       './BxB/methylKit/F78.Consensus.CpGs.methylKit.gz',
#                       './BxB/methylKit/F99.Consensus.CpGs.methylKit.gz')

ctrl.file.list <- list('./aligned2Angus/AxA/F53.Consensus.CpGs.methylKit.gz',
                        './aligned2Angus/AxA/F60.Consensus.CpGs.methylKit.gz',
                        './aligned2Angus/AxA/F7.Consensus.CpGs.methylKit.gz')

treat.file.list <- list('./aligned2Angus/BxB/F22.Consensus.CpGs.methylKit.gz',
                       './aligned2Angus/BxB/F65.Consensus.CpGs.methylKit.gz',
                       './aligned2Angus/BxB/F99.Consensus.CpGs.methylKit.gz')

file.list <- append(treat.file.list, ctrl.file.list)

############################### BP-DMR Analysis ################################
methObj <- methRead(file.list,
                    sample.id = list("F103","F105","F52","F53","F60","F7",
                                     "F22","F46","F56","F65","F78","F99"),
                    assembly='UOA-Brahman',
                    treatment=c(1,1,1,1,1,1,
                                0,0,0,0,0,0),
                    context='CpG',
                    resolution = 'base')#,
                    #dbtype="tabix",
                    #dbdir='./BBvsAA_base')

methObj <- methRead(file.list,
                    sample.id = list("F22","F65","F99",
                                     "F53","F60","F7"),
                    assembly='UOA-Angus',
                    treatment=c(1,1,1,
                                0,0,0),
                    context='CpG',
                    resolution = 'base')

meth <- filterByCoverage(methObj,
                                     lo.count = 10, lo.perc=NULL,
                                     hi.count = NULL, hi.perc = 99.9)#,
                                     #save.db = TRUE)

meth <- normalizeCoverage(meth)

meth <- unite(meth, destrand = T,
              min.per.group = 3L,
              mc.cores = 6)#,
              #save.db = T,
              #suffix = 'min5_destrand')

makeMethylDB(meth,'./aligned2Angus/AAvsBB.Consensus.FEMALE.base')

meth <- readMethylDB('./aligned2Angus/AAvsBB.Consensus.FEMALE.base/methylBase_min5_destrand.txt.bgz')
#meth <- readMethylDB('./BBvsAA.consensus.base/methylBase_min5_destrand.txt.bgz')

# Set male = 0, female = 1
#covariates <- data.frame(sex=c(0,0,0,1,1,1,
#                               1,0,0,1,0,1))

meth_Diff <- calculateDiffMeth(meth,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth_bp',
                               slim = T)

meth_tiles <- tileMethylCounts(meth,
                               mc.cores = 6,
                               save.db=TRUE,
                               suffix = 'tiled')

# Set male = 0, female = 1
covariates <- data.frame(sex=c(0,0,0,1,1,1,
                               1,0,0,1,0,1))

meth_Diff <- calculateDiffMeth(meth_tiles,
                               mc.cores = 6,
                               save.db = T,
                               suffix = 'diffMeth',
                               slim = T)

dmrs <- getMethylDiff(meth_Diff, type='hyper', difference=10, save.db=F)
dim(dmrs)


meth <- readMethylDB('./BBvsAA.base/methylDiff_min5_destrand_tiled_w1000s500_diffMeth.txt.bgz')

up_dmrs <- getMethylDiff(meth, difference = 10, type='hyper',save.db = F)
down_dmrs <- getMethylDiff(meth, difference = 10, type='hypo',save.db = F)

up_dmrs_granges <- as(up_dmrs, "GRanges")
up_dmrs_df <- data.frame(seqnames=seqnames(up_dmrs_granges),
                      starts=start(up_dmrs_granges)-1,
                      ends=end(up_dmrs_granges),
                      names=up_dmrs_granges$meth.diff,
                      scores=up_dmrs_granges$qvalue,
                      strands=strand(up_dmrs_granges))
options(scipen = 100, digits=4)
write.table(up_dmrs_df, file = './BBvsAA.base/BBvsAA_DMRs_10pHyper.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

down_dmrs_granges <- as(down_dmrs, "GRanges")
down_dmrs_df <- data.frame(seqnames=seqnames(down_dmrs_granges),
                         starts=start(down_dmrs_granges)-1,
                         ends=end(down_dmrs_granges),
                         names=down_dmrs_granges$meth.diff,
                         scores=down_dmrs_granges$qvalue,
                         strands=strand(down_dmrs_granges))
options(scipen = 100, digits=4)
write.table(up_dmrs_df, file = './BBvsAA.base/BBvsAA_DMRs_10pHypo.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

library(genomation)
gene_bed <- readTranscriptFeatures('../../Reference_Genomes/common_UOA-Brahman-Angus/Brahman_oriented2_ARS/Bos_indicus_hybrid.UOA_Brahman_1.104.FLIPPED.bed')
annotateWithGeneParts(as(up_dmrs, "GRanges"),gene_bed)

tss_df <-data.frame(seqnames=seqnames(gene_bed$TSSes),
           starts=start(gene_bed$TSSes)-1,
           ends=end(gene_bed$TSSes),
           names=gene_bed$TSSes$name,
           scores=gene_bed$TSSes$score,
           strands=strand(gene_bed$TSSes))

write.table(tss_df, file = '../../Reference_Genomes/common_UOA-Brahman-Angus/Brahman_oriented2_ARS/Bos_indicus_hybrid.UOA_Brahman_1.104.FLIPPED.TSSes.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

prom_df <-data.frame(seqnames=seqnames(gene_bed$promoters),
                    starts=start(gene_bed$promoters)-1,
                    ends=end(gene_bed$promoters),
                    names=gene_bed$promoters$name,
                    scores=gene_bed$promoters$score,
                    strands=strand(gene_bed$promoters))

write.table(prom_df, file = '../../Reference_Genomes/common_UOA-Brahman-Angus/Brahman_oriented2_ARS/Bos_indicus_hybrid.UOA_Brahman_1.104.FLIPPED.PROMOTERS.bed',
            quote = F, sep='\t',row.names = F, col.names = F)

tsses <- gene_bed$TSSes
prom <- gene_bed$promoters
