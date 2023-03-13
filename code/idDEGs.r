library(edgeR)

seqdata <- read.delim('./all.aligned2Brahman.RNA.txt', stringsAsFactors = FALSE) # Replace with aligned to Angus as required.
sampleinfo <- read.delim('./sampleinfo.txt', stringsAsFactors = FALSE)
sampleinfo$Genetics <- as.factor(sampleinfo$Genetics)
sampleinfo$Sex <- as.factor(sampleinfo$Sex)
sampleinfo$Batch <- as.factor(sampleinfo$Batch)


# Rename the columns so they look nicer
colnames(seqdata) <- gsub('.sorted.bam', '', colnames(seqdata))
colnames(seqdata) <- gsub('X', 'F', colnames(seqdata))

seqdata$F46 <- NULL
seqdata$F52 <- NULL

# Check if counts data and sample info are in the same order
table(colnames(seqdata) == sampleinfo$SampleName)

# Change order of counts data so that it matches the sample info
colorder <- c(as.character(sampleinfo$SampleName))
countsdata <- seqdata[, colorder]
table(colnames(countsdata) == sampleinfo$SampleName)




x <- DGEList(countsdata, samples=sampleinfo)
x
names(x)

sel <- rowSums(cpm(x$counts)>0.5)>=3#
x$counts <- x$counts[sel,]#

write.csv(cpm(x$counts), file='../methylation/methPipe/aligned2Brahman/Brahman.gene.expression.CPM.csv')
x

genetics <- as.factor(sampleinfo$Genetics)#
sex <- as.factor(sampleinfo$Sex)#
cols <- c(rep('red',6), rep('green',6), rep('orange', 6), rep('blue',6))

cols <- rep("blue",ncol(x$counts))
cols[sex==2] <- "red"

dim(x)

x <- calcNormFactors(x, method="TMM")#
batch <- sampleinfo$Batch
des <- model.matrix(~0+genetics+sampleinfo$Batch)
row.names(des) <- sampleinfo$SampleName
colnames(des) <- c('BIBI','BIBT','BTBI','BTBT','batch2')
des

x <- estimateDisp(x,des)
v <- voomWithQualityWeights(x,design=des,plot=F,col=as.numeric(genetics))

plotMDS(v,label=genetics,col=cols,dim.plot=c(1,2),main="")
v2=removeBatchEffect(v,batch=batch)#
legend("left",legend=c("BIBI","BIBT",'BTBT', 'BTBT'),text.col=c("red",'green','orange','blue'))
tiff('/Users/callummacphillamy/OneDrive/PhD/brahman_v_angus/MANUSCRIPT/figures/gene_expression_mds_Brahman.tiff',
     width=10,height=10,units='in',res = 400)
plotMDS(v,label=genetics,col=cols,dim.plot=c(1,2),main="")
legend("left",legend=c("BIBI","BIBT",'BTBT', 'BTBT'),text.col=c("red",'green','orange','blue'))
dev.off()


plotMDS(v2,label=genetics,col=cols,dim.plot=c(1,2),main="C-Transcriptome-Liver-Genetics")


# Here, if the logCPM is positive, then it represents 
# increased expression in left compared to right
# If the logCPM is negative, then it means it is more highly
# expressed in right compared to left.
contr=makeContrasts("BixBi-BixBt"=BixBi-BixBt,
                    "BixBi-BtxBi"=BixBi-BtxBi,
                    "BixBi-BtxBt"=BixBi-BtxBt,
                    "BixBt-BtxBi"=BixBt-BtxBi,
                    "BixBt-BtxBt"=BixBt-BtxBt,
                    "BtxBi-BtxBt"=BtxBi-BtxBt,
                        levels=des)

vfitgenetics <- lmFit(v)

colnames(vfitgenetics$coefficients) <- rownames(contr)
vfitcgenetics=contrasts.fit(vfitgenetics,contr)
vfitcgenetics = eBayes(vfitcgenetics)
resultsgenetics=decideTests(vfitcgenetics,p.value=0.05)
summary(decideTests(vfitcgenetics,p.value=0.05))

expression=v$E

comparisons <- colnames(vfitcgenetics$coefficients)

for (num in seq(1,6)) {
  print(num)
  print(comparisons[num])
  ttable <- topTable(vfitcgenetics,coef=num, n=Inf, p.value=0.05)
  write.csv(ttable, file=paste0('./DEGs/aligned2Brahman/Brahman.liver.',comparisons[num],'.pval_0.05.csv'))
  
}


