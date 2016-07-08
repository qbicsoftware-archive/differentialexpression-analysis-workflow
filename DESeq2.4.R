#Analysis pipeline from HTseq for differential expression analysis using DESeq2
# direct output to a file
con <- file("./logs/diff_expression_analysis_DESeq2.4.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

#load necessary libs
library("plyr")
library("ggplot2")
library("reshape2")
library("DESeq2")
library("gplots") 
library("RColorBrewer") 
library("genefilter") 
library("sqldf")
library("vsn")


#clean up
rm(list = ls(all = TRUE)) # clear all variables
graphics.off()

path <- getwd()
setwd(path)


######################################
# create directories needed
ifelse(!dir.exists("./result/DESeq2_run"), dir.create("./result/DESeq2_run"), FALSE)
dir.create("./result/DESeq2_run/raw_counts")
dir.create("./result/DESeq2_run/results")
dir.create("./result/DESeq2_run/results/plots")
dir.create("./result/DESeq2_run/results/tables")
dir.create("./result/DESeq2_run/results/final")
######################################

###1)check data
#put count table to counts/
files = list.files("./data/",pattern = "*Counts*")
for (i in files)
{
  cmd=paste("cp ", paste(i), "./result/DESeq2_run/raw_counts/",sep="" )
  system(cmd)
}

###2) Prepare count table
lof <- list.files("./result/DESeq2_run/raw_counts/")
lof <- lof[grep(".*.txt", lof)]  #extract only files ending on .txt
first.sample <- read.delim(paste0(lof[1]),header=F,row.names=1)
count.table <- data.frame(first.sample)   #writes an empty df before the for loops
names <- list(paste0(lof[1]))
for (s in lof [2:length(lof)]) {
   fname <- paste0(s)
   names <- c(names,fname)
   column <- read.delim(fname,header=F,row.names=1)
   count.table <- cbind(count.table,s=column) #binds the data frames together
}

#lof = list.files('./result/DESeq2_run/raw_counts',pattern = "counts*",full.names = T)
#count.table <- read.table(lof,  header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE,row.names=1)
#names(count.table)
colnames(count.table) <- names

#colnames(count.table) <- lof
#the fact that colnames are not informative can be disregarded 
#here as DESeq dataset will be built later anyhow....

###2.1) remove lines with "__" from HTSeq
count.table <- count.table[!grepl("__",row.names(count.table)),]
#dim(count.table)
#tail(count.table)


###3) load metadata
setwd(path)
metadata <- read.table("groups.txt", sep="\t", header=TRUE,na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE,row.names = 1)
sapply(metadata,class)
#make sure metadata is factor where needed
metadata$treatment <- as.factor(metadata$treatment)
metadata$time <- as.factor(metadata$time)
metadata$group <- as.factor(metadata$group)
metadata$bio_rep <- as.factor(metadata$bio_rep)

sapply(metadata,class)

#REMOVE a dataset as indicated during the QC pipeline
#metadata <- metadata[-5,]
#count.table <- count.table[,-5]

#check metadata and count table sorting
identical(names(count.table),row.names(metadata))
#write a check and stop function to stop script if that
#command does not lead to TRUE

###4) run DESeq function
#As group is the variable of interest, we put it at the end of the formula. Thus
cds <- DESeqDataSetFromMatrix( countData =count.table, colData =metadata, design =~group )
#cds <- DESeq(cds,parallel = TRUE)
cds <- DESeq(cds)
#cds
identical(colnames(cds),names(count.table))
#levels of interesting conditions to inspect:
groups <- levels(cds$group)
#treatment <- levels(cds$treatment)
#condition <- levels(cds$condition)

### 4.1) sizeFactors(cds) as indicator of library sequencing depth
sizeFactors(cds)
write.table(sizeFactors(cds),"./result/DESeq2_run/results/tables/sizeFactor_libraries.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = F, qmethod = c("escape", "double"))

#write raw counts to file
write.table(count.table, "./result/DESeq2_run/results/tables/raw.read.counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))
#the same as in assay(cds)....DE should be performed on raw read counts, but how to get the final values
#still problems to add ID for row.names

#if needed, also provide "normal" log2 counts, but add arbitray value  +0.1 to each value to avoid to take log2 from 0
log2counts <- log2(assay(cds)+0.1)
write.table(log2counts, "./result/DESeq2_run/results/tables/log2counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))

setwd(path)

#check some cook values to estimate cook distances between samples
assays(cds)[["cooks"]][1:10,]
#built in here a function that looks for Cook distances and 
#continues if there are no genes clearly falling out because 
#of being assigned to be outliers.
#Cooks distances: get important for example when checking knock-out and overexpression studies
pdf("./result/DESeq2_run/results/plots/Cooks-distances.pdf")
par(mar=c(10,3,3,3))
par( mfrow = c(1,2))
boxplot(log10(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log10-Cooks")
boxplot(log2(assays(cds)[["cooks"]]), range=0, las=2,ylim = c(-15, 15),main="log2-Cooks")
dev.off()


####4.2) contrasts
#to inspect
resultsNames(cds)

#contrasts to extract: (simple one way comparison)
normal_vs_tumor = results(cds, contrast=c("group","normal","tumor"))
head(normal_vs_tumor)


#make a list l and add results to it if more than one contrast is analyzed:
#l <- list(assign(paste("results", levels[1], levels[2], sep="_"), results(cds, contrast=c("genetic_background",levels[1],levels[2]),cooksCutoff=FALSE)))
l <- list(normal_vs_tumor=normal_vs_tumor)
length(l)


###4.3) significance after adjusted p values for multiple testing using BH
resSig <- sapply(l, function(i) FUN=i[ which(i$padj < 0.05), ])
#sort by FC
#to see the top positive DEGs:
head(sapply(resSig, function(i) FUN=i[ order(-i$log2FoldChange ), ]),n=5)
#to see the top negative DEGs:
head(sapply(resSig, function(i) FUN=i[ order(i$log2FoldChange ), ]),n=5)


#write l to file containing all unique combinations contrasts
df <- as.data.frame(l)
#quick summary
#sum(l$padj < 0.05, na.rm=TRUE ) 
#empty
head(df)
names(df)
write.table(df, "./result/DESeq2_run/results/final/final_list_DESeq2.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))


###Data transformation
#The function rlogTransform returns a SummarizedExperiment object which contains the rlog-transformed values in its assay slot:
rld <- rlog(cds)

##VST Transformation as done in DESeq:   
#according to DESeq package more sensitive to oultliers thus e.g. when library sizes vary widely. Then it would be better to trust the rlog data.
vsd <- varianceStabilizingTransformation(cds)

#example to compare rlog and VST
assay(rld)[c(1:3),]
assay(vsd)[c(1:3),]

write.table(assay(rld), "./result/DESeq2_run/results/tables/rlog_transformed.read.counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = NA, qmethod = c("escape", "double"))
write.table(assay(vsd), "./result/DESeq2_run/results/tables/vst_transformed.read.counts.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = NA, qmethod = c("escape", "double"))

#Use these values when you want to cluster samples or genes (trees, PCA plots, heatmaps) as they allow you to group samples and genes unbiases by technical artefacts such as sequencing depth.
#If your downstream method is designed to take as input expression matrices similar to normalized microarray datasets (log scale) then
#you can use rlog or VST, and use the matrix accessed with assay(object).For example the following extracted table would be good for GSEA

#################
##Diagnostic plots
setwd(path)
#The function plotDispEsts visualizes DESeqs dispersion estimates: 
pdf("./result/DESeq2_run/results/plots/Dispersion_plot.pdf")
plotDispEsts(cds, ylim = c(1e-5, 1e8))
dev.off()
png("./result/DESeq2_run/results/plots/Dispersion_plot.png")
plotDispEsts(cds, ylim = c(1e-5, 1e8))
dev.off()
#Effects of transformations on the variance
notAllZero <- (rowSums(counts(cds))>0) 
pdf("./result/DESeq2_run/results/plots/Effects_of_transformations_on_the_variance.pdf")
par(oma=c(3,3,3,3))
par(mfrow = c(1, 3))
meanSdPlot(log2(counts(cds,normalized=TRUE)[notAllZero,] + 1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()


#index the samples
normal = which(cds[["group"]] == "normal")
tumor = which(cds[["group"]] == "tumor")
#Scatter plot -- effect of trafo
pdf("./result/DESeq2_run/results/plots/scatter_plot_effect_of_transformation.pdf")
par( mfrow = c(1,2)) 
plot( log2( 1+counts(cds, normalized=TRUE)[, normal[1]:tumor[1]]), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, normal[1]:tumor[1]], col="#00000020", pch=20, cex=0.3 )
dev.off()
#Note that, in order to make it easier to see where several points 
#are plotted on top of each other, we set the plotting color to a 
#semi-transparent black (encoded as #00000020)and changed the points to 
#solid disks (pch=20)with reduced size (cex=0.3).


# scatter plot of rlog or vst transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld),n=2)
head(assay(vsd),n=2)

###Sample distances
#A useful first step in an RNA-Seq analysis is often to assess overall similarity between samples: Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment???s design?
#We use the R function dist to calculate the Euclidean distance between samples. To avoid that the distance measure is dominated by a few highly variable genes, and have a roughly equal contribution from all genes, we use it on the rlog-transformed data:

sampleDists <- dist(t(assay(rld)))  #Note the use of the function t to transpose the data matrix. We need this because dist calculates distances between data rows and our samples constitute the columns.
sampleDists #high values=high distance between samples, 0=no distance,darkblue here....
#We visualize the distances in a heatmap, using the function heatmap.2 from the gplots package.
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(colnames(rld),rld$genetic_background, sep="--")
#colnames(sampleDistMatrix) <- paste(colnames(rld),rld$genetic_background, sep="--")
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255) 

###1) Visualization of distance using Heatmaps
#Heatmap of Euclidean sample distances after rlog transformation
pdf("./result/DESeq2_run/results/plots/Heatmaps_of_distances.pdf")
par(oma=c(3,3,3,3))
heatmap.2(sampleDistMatrix, trace="none", col=colours,cexRow=0.5,cexCol=0.5,keysize=1.5,key.title=NA,key.xlab="value",key.ylab="count",key.par=list(cex=0.6))
dev.off()
png("./result/DESeq2_run/results/plots/Heatmaps_of_distances.png")
par(oma=c(3,3,3,3))
heatmap.2(sampleDistMatrix, trace="none", col=colours,cexRow=0.5,cexCol=0.5,keysize=1.5,key.title=NA,key.xlab="value",key.ylab="count",key.par=list(cex=0.6))
dev.off()

####2) Visualization of distance using PCA plots
#Another way to visualize sample-to-sample distances is a principal-components analysis (PCA). In this ordination method, the data points (i.e., here, the samples) are projected onto the 2D plane such that they spread out optimally. 
pdf("./result/DESeq2_run/results/plots/PCA_plot_of_distances.pdf")
print(plotPCA(rld,intgroup=c("group")))
dev.off()
#for html:
png("./result/DESeq2_run/results/plots/PCA_plot_of_distances.png")
print(plotPCA(rld,intgroup=c("group")))
dev.off()


##
###Gene clustering
# In the heatmap,the dendrogram at the side shows us a hierarchical clustering of the samples. Such a clustering can also be performed for the genes.  
# Since the clustering is only relevant for genes that actually carry signal, one usually carries it out only for a subset of most highly variable genes. Here, for demonstration, let us select the 35 genes with the highest variance across samples:
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 50)
#The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene???s average across all samples. Hence, we center and scale each genes??? values across samples, and plot a heatmap
pdf("./result/DESeq2_run/results/plots/heatmap_of_gene_clusters.pdf")
par(oma=c(3,3,3,3))
heatmap.2(assay(rld)[topVarGenes, ],scale="row",trace="none",dendrogram="row",col=colorRampPalette( rev(brewer.pal(9, "RdBu")))(255),cexRow=0.5,cexCol=0.5)
#take out labCol=paste(colnames(rld),rld$genetic_background, sep="--")
dev.off()
## trace="none"  turns off trace lines inside the heat map
#We can now see blocks of genes which covary across samples.
##

#better than to show genes with highest variance see above go back to the specific gene lists from resSig...did not do it here yet

# #select ctnnb1 vs WT_T which were also not found by comparing ctnnb1 vs cl or any other contrast
# # e.g. do
# assay(rld)[ctnnb1_vs_WT_T_specific[,1],]
# pdf("results/plots/heatmap_of_ctnnb1_vs_WT_T.pdf")
# par(oma=c(3,3,3,3))
# heatmap.2(assay(rld)[ctnnb1_vs_WT_T_specific[,1],],scale="row",trace="none",dendrogram="column",col=colorRampPalette( rev(brewer.pal(9, "RdBu")))(255),labCol=paste(colnames(rld),rld$condition, sep="--"),cexRow=0.7,cexCol=0.7)
# dev.off()

#MA plots to visualize DE
#plot "DESeqResults" using plotMA function. For that results(cds) should be written again as DESeqResults:
pdf("./result/DESeq2_run/results/plots/all_results_MA_plot.pdf")
plotMA(results(cds),ylim = c(-4, 4))
dev.off()
#plot to html file
png("./result/DESeq2_run/results/plots/all_results_MA_plot.png")
plotMA(results(cds),ylim = c(-4, 4))
dev.off()

#MA-plots of tests of log2 fold change vs mean expression value to show dependency of calls to be DE with read count numbers.  
#the gene at the higher end (log2fold change >5= **ERF6_1**, the one at the lower end of the plot was VIT_03s0088g00690, which is annotated with Jerome's paper as a Pathogenesis-related.protein.1.precursor.(PRP.1).

#The MA plot above highlights an important property of RNA-Seq data. For weakly expressed genes, we have no chance of seeing differential expression, because the low read counts suffer from so high Poisson noise that any biological effect is drowned in the uncertainties from the read counting. We can also show this by examining the ratio of small p values (say, less than, 0.01) for genes binned by mean normalized count:
#  to check dependency of small p values from mean normalized counts: rationale to do so: for weakly expressed genes, there is no chance of seeing differential expression, because the low read counts suffer from so high Poisson noise that any biological effect is drowned in the uncertainties from the read counting.

setwd(path)


#Multiple hypothesis testing
# 1) create bins using the quantile function 
qs <- c( 0, quantile(results(cds)$baseMean[results(cds)$baseMean > 0], 0:4/4 )) #The generic function quantile produces sample quantiles corresponding to the given probabilities. The smallest observation corresponds to a probability of 0 and the largest to a probability of 1.
# 2) "cut" the genes into the bins 
bins <- cut(results(cds)$baseMean, qs ) 
# rename the levels of the bins using the middle point 
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)])) 
# 3) calculate the ratio of p values less than .01 for each bin 
ratios <- tapply(results(cds)$pvalue, bins, function(p) mean(p < .01, na.rm=TRUE )) 
# 4) plot these ratios 
pdf("./result/DESeq2_run/results/plots/dependency_small.pval_mean_normal.counts.pdf")
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
dev.off()

#At first sight, there may seem to be little benefit in filtering out these genes with low mean normalized count. After all, the test found them to be non-significant anyway. However, these genes have an influence on the multiple testing adjustment, whose performance improves if such genes are removed. By removing the weakly-expressed genes from the input to the FDR procedure, we can find more genes to be significant among those which we keep, and so improved the power of our test. This approach is known as independent filtering.


# pdf("DESeq2/results/plots/number.of.rejections.pdf")
# plot(attr(results(cds),"filterNumRej"),type="b", xlab="quantiles of 'baseMean'=mean of normallized counts (theta)", ylab="number of rejections")
# dev.off()

#pdf("DESeq2/results/plots/number.of.rejections.pdf")
#plot(metadata(results(cds))$filterNumRej,
#     type="b", ylab="number of rejections",
#     xlab="quantiles of filter")
#lines(metadata(results(cds))$lo.fit, col="red")
#abline(v=metadata(results(cds))$filterTheta)
#dev.off()

##This plot shows the process of Independent filtering. DESeq2 automatically determines a threshold, filtering on mean normalized count, which maximizes the number of genes which will have an adjusted p value less than a critical value. Here we will get ~XXX genes with an adjusted p value smaller than a critical value (theta)
#The term independent highlights an important caveat. Such filtering is permissible onl if the filter criterion is independent of the actual test statistic. Otherwise, the filtering would invalidate the test and consequently the assumptions of the BH procedure. This is why we filtered on the average over all samples: this filter is blind to the assignment of samples to the treatment and control group and hence independent.

#Histogram of passed and rejected hypothesis for contrast ctnnb1 vs control liver
use <- results(cds)$baseMean > metadata(results(cds))$filterThreshold
table(use)
h1 <- hist(results(cds)$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(results(cds)$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c('do not pass'="khaki", 'pass'="powderblue")
pdf("./result/DESeq2_run/results/plots/histogram_of_p.values.pdf")
barplot(height = rbind(h1$density, h2$density), beside = FALSE,
        col = colori, space = 0, main = "", xlab="p value",ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topleft", fill=rev(colori), legend=rev(names(colori)))
dev.off()

#The p value histogramm here shows that filtering did not take out visibly many genes 
#but when looking at the table(use)
table(use)
#we see that XXX genes failed to pass independent filtering. This will ameliorate the multiple 
#testing problemand thus the severity of a multiple testing adjustment 
#by removing a background set of hypotheses whose p values are not equally distributed here,
#but around 1. It does not make sense to arbitrarily increase the filter threshold by e.g. 
#saying 
#use <- results$baseMean > 1 to indicate that baseMean should be bigger than 1 
#as this will not improve number of genes being detected after multiple testing using FDR 
#(see plot "number.of.rejections.pdf" when walking to the right...)

sessionInfo()

#END
#FINISHED HERE

