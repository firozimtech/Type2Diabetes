if(FALSE){
#########################################################
#  Getting the "Microarray Raw" data directly from GEO  #
#########################################################
library(GEOquery)
getGEOSuppFiles("GSE2899")


#########################################################
#  Create directory "data" under your working directory #
#  Unpack CEL file  from "GSE20986_RAW.tar" to "data"   #  
#  Decompressed the CEL file
#########################################################
# NOTE: R on OS X 10.9 (Mavericks) seems to set a wrong TAR environment variable,you can fix this by executing:
Sys.setenv(TAR = '/usr/bin/tar')

untar("GSE20986/GSE20986_RAW.tar", exdir="data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
cels

#########################################################
#  Create a "phenodata.txt", 3 tab-delimited columns    #
#  named ‘Name’, ‘FileName’ and ‘Target’.               #
#########################################################

}

#########################################################
#       Loading and Normalizing the data		#
#########################################################
library(simpleaffy)
celfiles <- read.affy(covdesc="phenodata.txt", path="data")
celfiles
celfiles.gcrma <- gcrma(celfiles)
# NOTE: Expression values by RMA or GCRMA are  log2 transformed. Amount of mRNA is roughly 2x, where x is the expression value
# signal intensity values generated from MAS5.0 is not log2 transformed

#########################################################
#       Quality control Checks using probe intensity    #
#          before and after normalization               #
#########################################################
library(RColorBrewer)
cols <- rainbow(12)

# Boxplot of unnormalized intensity values

pdf("UnNormalizedGraph.pdf",width=7,height=5)
# If plot is not fit into window
par(mar = c(9,4,4,2) + 0.1)
boxplot(celfiles, col=cols, las=2)
dev.off()

# Boxplot of normalized intensity values
library(affyPLM)
pdf("NormalizedGraph.pdf",width=7,height=5)
# If plot is not fit into window
par(mar = c(9,4,4,2) + 0.1)
boxplot(celfiles.gcrma, col=cols, las=2)
dev.off()


# Plot a density vs. log intensity histogram for the unnormalized data
pdf("UnNormalizedHistogram.pdf",width=7,height=5)
# If plot is not fit into window
par(mar = c(9,4,4,2) + 0.1)
hist(celfiles, col=cols)
dev.off()

# Plot a density vs. log intensity histogram for the normalized data
pdf("NormalizedHistogram.pdf",width=7,height=5)
# If plot is not fit into window
par(mar = c(9,4,4,2) + 0.1)
hist(celfiles.gcrma, col=cols)
dev.off()

# Relationships between the samples(Normalized) using hierarchical clustering:
eset <- exprs(celfiles.gcrma)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
pdf("NormalizedSampleCluster.pdf",width=7,height=5)
# If plot is not fit into window
par(mar = c(9,4,4,2) + 0.1)
plot(clusters)
dev.off()


###############################################################
# Filter out uninformative data (a) Control probesets         #
# (b) internal controls, (c) removing genes with low variance,#
# which will be unlikely to pass statistical tests for        #
# differential expression, or are expressed uniformly close to#
# background detection levels.                                #
###############################################################

#celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=TRUE,remove.dupEntrez=TRUE)
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE,remove.dupEntrez=FALSE)
# What got removed and why?
celfiles.filtered$filter.log

filterEset <- exprs(celfiles.filtered$eset)
dim(filterEset)

#########################################################################################
# Principal component analysis (PCA): dimension reduction and data visualization method #
# The axes in the PCA map represent 3 major components identified. It shows us the      #
# overall gene expression patterns. Here we use PCA to visualize samples relationship   #
# and to detect outliers.                                                               #
#########################################################################################

library(rgl)
pca <- prcomp(t(filterEset), scale=TRUE)
summary(pca)

# Save summary of PCA 
sink("PCA_summary.txt")
print(summary(pca))
sink()

# Save all PCA without maximum row limit
options(max.print = 999999999)
sink("PCA_all.txt")
print(pca)
sink()


# choroid is red, huvec is green, iris is blue, retina is yellow
# Put the color in order of "name.CEL" file in phenodata.txt
#myColors <- c("Blue", "yellow", "yellow", "Blue", "yellow", "Blue",rep("Red", 3), rep("Green", 3))


myColors <- c("violet", "violet", "red", "red", "pink", "pink","green", "green", "gray", "gray","blue", "blue", "orange", "orange","purple", "purple")   

#myColors <- c("red", "red", "pink", "pink", "green", "green", "aquamarine", "aquamarine", "blue", "blue", "skyblue","skyblue", "yellow", "yellow", "darkgoldenrod", "darkgoldenrod")

plot3d(pca$x[, 1:3], col=myColors, xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")

legend3d("topright", legend = paste(c('Control_Adipose','Control_Adipose', 'T2D_Adipose', 'T2D_Adipose', 'Control_Liver','Control_Liver', 'T2D_Liver', 'T2D_Liver', 'Control_Muscle','Control_Muscle', 'T2D_Muscle', 'T2D_Muscle', 'Control_Islet','Control_Islet', 'T2D_Islet', 'T2D_Islet')), pch = 16, col = myColors, cex=2, inset=c(0.02))
# Save the above image into a file
rgl.postscript("PCA.pdf", fmt="pdf", drawText=TRUE)



##################### PCA for Adipose ##################
filterEsetAdipose <- filterEset[,1:4]
pcaAdipose <- prcomp(t(filterEsetAdipose), scale=TRUE)
summary(pcaAdipose)

# Save summary of PCA
sink("PCA_summary_adipose.txt")
print(summary(pcaAdipose))
sink()

# Save all PCA without maximum row limit
options(max.print = 999999999)
sink("PCA_adipose.txt")
print(pcaAdipose)
sink()

myColorsAdipose <- c("violet", "violet", "red", "red")

#myColorsAdiposeGrp<-factor(myColorsAdipose)

plot3d(pcaAdipose$x[, 1:3], col=myColorsAdipose, xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
legend3d("topright", legend = paste(c('Control_Adipose','Control_Adipose', 'T2D_Adipose', 'T2D_Adipose')), pch = 16, col = myColorsAdipose, cex=2, inset=c(0.02))

# Save the above image into a file
rgl.postscript("PCA_Adipose.pdf", fmt="pdf", drawText=TRUE)

#plot3d(pcaAdipose$x[, 1:3], col=as.integer(myColorsAdiposeGrp), xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
#legend3d("topright", legend = paste(c('Control_Adipose', 'T2D_Adipose')), pch = 16, col = myColorsAdipose, cex=2, inset=c(0.02))

##################### PCA for Liver ##################
filterEsetLiver <- filterEset[,5:8]
pcaLiver <- prcomp(t(filterEsetLiver), scale=TRUE)
summary(pcaLiver)

# Save summary of PCA
sink("PCA_summary_Liver.txt")
print(summary(pcaLiver))
sink()

# Save all PCA without maximum row limit
options(max.print = 999999999)
sink("PCA_Liver.txt")
print(pcaLiver)
sink()

myColorsLiver <- c("pink", "pink","green", "green")
plot3d(pcaLiver$x[, 1:3], col=myColorsLiver, xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
legend3d("topright", legend = paste(c('Control_Liver','Control_Liver', 'T2D_Liver', 'T2D_Liver')), pch = 16, col = myColorsLiver, cex=2, inset=c(0.02))

# Save the above image into a file
rgl.postscript("PCA_Liver.pdf", fmt="pdf", drawText=TRUE)




##################### PCA for Muscle ##################
filterEsetMuscle <- filterEset[,9:12]
#pcaMuscle <- prcomp(t(filterEsetMuscle), scale=TRUE)
#summary(pcaMuscle)

#### I m getting following error after above command #############
# Error in prcomp.default(t(filterEsetMuscle), scale = TRUE) : cannot rescale a constant/zero column to unit variance
# It says that you cannot rescale a column with zero variance to unit variance
# To convert to unit variance, you divide the values in the column by the standard deviation of the column. If the column has variance zero, the standard deviation is zero
# This also means you have genes (rows) that don't change expression across samples (identical values across samples).
# Anyway, a gene with no variability is de facto uninteresting, and removing those genes is a good idea.

#### To fix above problem ################ 

library(genefilter)
filterEsetMuscle <- filterEset[,9:12]
ind <- rowVars(filterEsetMuscle) < .Machine$double.eps
filterEsetMuscle<-filterEsetMuscle[!ind,]

# The process will remove three proble id from sample (100381_at, 161361_s_at, 99667_at)
pcaMuscle <- prcomp(t(filterEsetMuscle), scale=TRUE)
summary(pcaMuscle)



# Save summary of PCA
sink("PCA_summary_Muscle.txt")
print(summary(pcaMuscle))
sink()

# Save all PCA without maximum row limit
options(max.print = 999999999)
sink("PCA_Muscle.txt")
print(pcaMuscle)
sink()

myColorsMuscle <- c("gray", "gray","blue", "blue")
plot3d(pcaMuscle$x[, 1:3], col=myColorsMuscle, xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
legend3d("topright", legend = paste(c('Control_Muscle','Control_Muscle', 'T2D_Muscle', 'T2D_Muscle')), pch = 16, col = myColorsMuscle, cex=2, inset=c(0.02))

# Save the above image into a file
rgl.postscript("PCA_Muscle.pdf", fmt="pdf", drawText=TRUE)


##################### PCA for Islet ##################
filterEsetIslet <- filterEset[,13:16]
pcaIslet <- prcomp(t(filterEsetIslet), scale=TRUE)
summary(pcaIslet)

# Save summary of PCA
sink("PCA_summary_Islet.txt")
print(summary(pcaIslet))
sink()

# Save all PCA without maximum row limit
options(max.print = 999999999)
sink("PCA_Islet.txt")
print(pcaIslet)
sink()

myColorsIslet <- c("orange", "orange","purple", "purple")
plot3d(pcaIslet$x[, 1:3], col=myColorsIslet, xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
legend3d("topright", legend = paste(c('Control_Islet','Control_Islet', 'T2D_Islet', 'T2D_Islet')), pch = 16, col = myColorsIslet, cex=2, inset=c(0.02))

# Save the above image into a file
rgl.postscript("PCA_Islet.pdf", fmt="pdf", drawText=TRUE)

legend3d("topright", legend = paste(c(,'Control_Islet', 'T2D_Islet',)), pch = 16, col = myColorsIslet, cex=1, inset=c(0.02))




############ two D PCA with labeling ###########
#groups <- factor(rownames(pcaAdipose$x))
#plot(pcaAdipose$x[,1:2], col=groups)
#legend(0,0,groups, col=groups, pch=1)
#text(pcaAdipose$x[,1:2], labels=groups, pos=2)
############ END ##################



#####################################################
# Finding differentially expressed probesets:       #
# (A) extract information about the samples.        #
#####################################################
samples <- celfiles.gcrma$Target
samples

# convert into factors
samples <- as.factor(samples)
samples

# set up the experimental design

design <- model.matrix(~0 + samples)

## NOTE: Change header name according to your experiment phenodata
# colnames(design) <- c("choroid", "huvec", "iris", "retina")
colnames(design) <- c("BTRB_Adipose", "BTBR_Islet", "BTRB_Liver", "BTBR_Muscle", "J_Adipose", "J_Islet", "J_Liver", "J_Muscle")


# inspect the experiment design
design

#################################################################################
# Now we have normalized, filtered the data, and added a description of the     #
# data and experimental design. This will be consumed by the limma packages for #
# the DEG analysis.                                                             #
#################################################################################


#####################################################
# Finding differentially expressed probesets:       #
# (A) DEG analysis.                                 #
#####################################################

library(limma)
# fit the linear model to the filtered expression set
fit <- lmFit(filterEset, design)
## NOTE: Change following according to your experiment
#contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid,huvec_retina = huvec - retina, huvec_iris <- huvec - iris, levels=design)

contrast.matrix <- makeContrasts(BTRB_Adipose - J_Adipose, BTRB_Liver - J_Liver, BTBR_Muscle - J_Muscle, BTBR_Islet - J_Islet, levels=design )
# Now the contrast matrix is combined with the per-probeset linear model fit
huvec_fits <- contrasts.fit(fit, contrast.matrix)


############################################################################
# Calculating the differential expression by empirical Bayes shrinkage     #
# of the standard errors towards a common value, by computing              #
# (a) moderated t-statistics, (b) moderated F-statistic, (3) and log-odds: #
############################################################################

huvec_ebFit <- eBayes(huvec_fits)

# return the top 10 genes(probes)  for any given contrast
# coef=1 is huvec_choroid, coef=2 is huvec_retina, coef=3 is huvec_iris

topTable(huvec_ebFit, number=10, coef=1)

# Probes having absolute value of fold changes > 5
topTable(huvec_ebFit, coef=1, number=10000, lfc=5)

# To count the number of probes
nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=5))


# List with differentially expressed probes filtered by fold change > 5 and p-value < 0.05, huvec vs. choroid
probeset.list <- topTable(huvec_ebFit, coef=1, p.value=0.05, lfc=5)



#############################
# Tables of result files    #
#############################

# expression level (log2)of all cel files with all proble id (54675 affyids)[normalized]
write.table(eset, "NonfilterExprLog2Cel.gct", sep="\t", quote=FALSE)
nrows(eset)

# Table of expression log2 value of every samples after filter
write.table(filterEset, "ExprLog2Cel.txt", sep="\t", quote=FALSE)

# Table of expression log Fold Change value of every samples
write.table(huvec_ebFit, "AllResult.txt", sep="\t", quote=FALSE)

# Complete list of genes with p-values and fold change
# Coef=1, so we are just looking at "BTRB_Adipose - J_Adipose"
gene_list_1 <- topTable(huvec_ebFit, coef=1, number=1000000, sort.by="logFC")
write.table(gene_list_1, "resultsLogFC_BTRBAdipose-JAdipose.txt", sep="\t", quote=FALSE)

gene_list_2 <- topTable(huvec_ebFit, coef=2, number=1000000, sort.by="logFC")
write.table(gene_list_2, "resultsLogFC_BTRBLiver-JLiver.txt", sep="\t", quote=FALSE)

gene_list_3 <- topTable(huvec_ebFit, coef=3, number=1000000, sort.by="logFC")
write.table(gene_list_3, "resultsLogFC_BTRBMuscle-JMuscle.txt", sep="\t", quote=FALSE)

gene_list_4 <- topTable(huvec_ebFit, coef=4, number=1000000, sort.by="logFC")
write.table(gene_list_4, "resultsLogFC_BTRBIslet-JIslet.txt", sep="\t", quote=FALSE)



############################################################################
#      Annotating probe id with associated gene symbols                    #
############################################################################
#library(hgu133plus2.db)
library(mgu74av2.db)
library(annotate)
# To list the kinds of things that can be used as keys
# use the keytypes method
keytypes(mgu74av2.db)


# find and extract the SYMBOL, ENTREZID, GO ids associated with the  PROBEID
# GenOnt<-select(mgu74av2.db, rownames(gene_list_1), c("SYMBOL", "ENTREZID", "GO"), "PROBEID")
# use GO.db to find the Terms associated with those GOIDs
# library("GO.db")
# GenOnt2<-select(GO.db, GenOnt$GO, "TERM", "GOID")


# find and extract the SYMBOL, ENTREZID ids associated with the  PROBEID
GenAnn_1<-select(mgu74av2.db, rownames(gene_list_1), c("SYMBOL", "ENTREZID"), "PROBEID")
#res_1 <- cbind(gene_list_1,GenAnn_1)
write.table(GenAnn_1, "GenAnn_BTRBAdipose-JAdipose.txt", sep="\t", quote=FALSE)


GenAnn_2<-select(mgu74av2.db, rownames(gene_list_2), c("SYMBOL", "ENTREZID"), "PROBEID")
#res_2 <- cbind(gene_list_2,GenAnn_2)
write.table(GenAnn_2, "GenAnn_BTRBALiver-JLiver.txt", sep="\t", quote=FALSE)

GenAnn_3<-select(mgu74av2.db, rownames(gene_list_3), c("SYMBOL", "ENTREZID"), "PROBEID")
#res_3 <- cbind(gene_list_3,GenAnn_3)
write.table(GenAnn_3, "GenAnn_BTRBMuscle-JMuscle.txt", sep="\t", quote=FALSE)

GenAnn_4<-select(mgu74av2.db, rownames(gene_list_4), c("SYMBOL", "ENTREZID"), "PROBEID")
#res_4 <- cbind(gene_list_4,GenAnn_4)
write.table(GenAnn_4, "GenAnn_BTRBIslet-JIslet.txt", sep="\t", quote=FALSE)


#gene.symbols <- getSYMBOL(rownames(probeset.list), "hgu133plus2")
#results <- cbind(probeset.list, gene.symbols)
#write.table(results, "AllResults_Annotated.txt", sep="\t", quote=FALSE)


# The head command prints the first few values of a vector
head(gene_list_1$logFC)
head(gene_list_1$P.Value)
head(gene_list_1$adj.P.Val)


# Number of genes that have an absolute fold change greater than 2 and a p-value less than 0.05
sum(abs(gene_list_1$logFC) > 2 & gene_list_1$P.Value < 0.05)

####################################################################
# Volcano plot
###################################################################

require(ggplot2)

# Cut-off based upon Bonferroni
# Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
# Bonferroni cut-off is approximately 0.05/#tests_gene.
# gene_list_1$threshold = as.factor(abs(gene_list_1$logFC) > 2 & gene_list_1$P.Value < 0.05/no_of_genes)
# gene_list_1$threshold = as.factor(abs(gene_list_1$logFC) > 2 & gene_list_1$P.Value < 0.05/sum(abs(gene_list_1$logFC) > 2 & gene_list_1$P.Value < 0.05))
# gene_list_2$threshold = as.factor(abs(gene_list_2$logFC) > 2 & gene_list_2$adj.P.Val < 0.05)

# Cut-off based upon FDR
##Highlight genes that have an absolute fold change >= 2 and FDR< 0.05
# gene_list_1$threshold = as.factor(abs(gene_list_1$logFC) >= 2 & gene_list_1$adj.P.Val < 0.05)

##Construct the plot object
# g = ggplot(data=gene_list_1, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
#  geom_point(alpha=0.4, size=1.75) +
#  theme(legend.position = "none") +
#  xlim(c(-10, 10)) + ylim(c(0, 15)) +
#  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("Adipose") + theme(plot.title = element_text(hjust = 0.5))
# g
# ggsave("VolcanoBTRBAdipose-JAdipose.pdf", width=7,height=5)


# Cut-off based upon P-value

gene_list_1$threshold = as.factor(abs(gene_list_1$logFC) >= 2 & gene_list_1$P.Value < 0.05)

##Construct the plot object
g = ggplot(data=gene_list_1, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("Adipose") + theme(plot.title = element_text(hjust = 0.5))
g
ggsave("VolcanoBTRBAdipose-JAdipose.pdf", width=7,height=5)


gene_list_2$threshold = as.factor(abs(gene_list_2$logFC) >= 2 & gene_list_2$P.Value < 0.05)
##Construct the plot object
g = ggplot(data=gene_list_2, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("Liver") + theme(plot.title = element_text(hjust = 0.5))
g
ggsave("VolcanoBTRBLiver-JLiver.pdf", width=7,height=5)


gene_list_3$threshold = as.factor(abs(gene_list_3$logFC) >= 2 & gene_list_3$P.Value < 0.05)
##Construct the plot object
g = ggplot(data=gene_list_3, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("Muscle") + theme(plot.title = element_text(hjust = 0.5))
g
ggsave("VolcanoBTRBMuscle-JMuscle.pdf", width=7,height=5)



gene_list_4$threshold = as.factor(abs(gene_list_4$logFC) >= 2 & gene_list_4$P.Value < 0.05)
##Construct the plot object
g = ggplot(data=gene_list_4, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+ ggtitle("Islet") + theme(plot.title = element_text(hjust = 0.5))
g
ggsave("VolcanoBTRBIslet-JIslet.pdf", width=7,height=5)


##########
## find and extract the SYMBOL, ENTREX, GO ids associated with the  PROBEID
GenOn<-select(mgu74av2.db, rownames(gene_list_1), c("SYMBOL", "ENTREZID", "GO"), "PROBEID")
