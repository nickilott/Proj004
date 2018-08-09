###########################################################
###########################################################
###########################################################
# compare RNA-seq derived fold changes to array fold
# changes - for CPZ vs control comparison
###########################################################
###########################################################
###########################################################

library(ggplot2)
library(plyr)

readData <- function(x){

	 x <- read.csv(x, header=T, stringsAsFactors=F, sep="\t")
	 return(x)
	 }

# read data
array <- readData("CPZUninfected-UntreatedUnInfected.result")
rnaseq <- readData("cpz_vs_control.result")


# get unique genes on array - take the maximum fold change
d <- ddply(array, .(gene), summarize, maxf=max(logFC))

# get the intersection of gene names
common.genes <- intersect(d$gene, rnaseq$gene_name)

# subset the data
rnaseq <- rnaseq[rnaseq$gene_name %in% common.genes,]
d <- d[d$gene %in% common.genes,]

# merge
dat <- merge(rnaseq, d, by.x="gene_name", by.y="gene", all.x=T, all.y=T)
dat$sig <- ifelse(dat$padj < 0.05 & !(is.na(dat$padj)), "yes", NA)

# plot
cor1 <- cor(dat$log2FoldChange, dat$maxf, use="pairwise.complete.obs")

p1 <- ggplot(dat, aes(x=log2FoldChange, y=maxf))
p2 <- p1 + geom_point(pch=18, colour="grey")
p3 <- p2 + theme_bw()
p4 <- p3 + xlim(c(-5,5)) + ylim(c(-5,5))
p5 <- p4 + annotate(geom="text", -5, 5, label=paste("r=", round(cor1,2), sep=""))
p6 <- p5 + geom_smooth(data=dat, aes(x=log2FoldChange, y=maxf), method="lm", inherit.aes=F)
p7 <- p6 + geom_point(dat=dat, aes(x=log2FoldChange, y=maxf, colour=sig), pch=18) + scale_colour_manual(values=c("red4"))
ggsave("scatter_fold_cpz.pdf")

###########################################################
###########################################################
###########################################################
# do the same in salmonella infection
###########################################################
###########################################################
###########################################################

# read data
array <- readData("UntreatedInfected-UntreatedUnInfected.result")
rnaseq <- readData("sal_vs_control.result")

# get sig.genes in array
sig.array <- unique(array$gene[array$adj.P.Val < 0.05])
not.sig.array <- setdiff(array$gene, sig.array)

# get unique genes on array - take the maximum fold change
d <- ddply(array, .(gene), summarize, maxf=max(logFC))

# get the intersection of gene names
common.genes <- intersect(d$gene, rnaseq$gene_name)

# subset the data
rnaseq <- rnaseq[rnaseq$gene_name %in% common.genes,]
d <- d[d$gene %in% common.genes,]

# merge
dat <- merge(rnaseq, d, by.x="gene_name", by.y="gene", all.x=T, all.y=T)

dat$label <- ifelse(dat$padj < 0.05 & !(is.na(dat$padj)) & dat$gene_name %in% sig.array, "Both", "Neither")
dat$label <- ifelse(dat$padj < 0.05 & !(is.na(dat$padj)) & dat$gene_name %in% not.sig.array, "RNA-seq", dat$label)
dat$label <- ifelse(dat$padj > 0.05 & !(is.na(dat$padj)) & dat$gene_name %in% sig.array, "Microarray", dat$label)
dat$label <- ifelse(dat$padj > 0.05 & !(is.na(dat$padj)) & dat$gene_name %in% not.sig.array, "Neither", dat$label)
dat$label <- ifelse(is.na(dat$label), "Neither", dat$label)
dat$label <- factor(dat$label, levels=c("Neither", "Both", "Microarray", "RNA-seq"))

nboth <- nrow(dat[dat$label == "Both",])
nneither <- nrow(dat[dat$label == "Neither",])
narray <- nrow(dat[dat$label == "Microarray",])
nrnaseq <- nrow(dat[dat$label == "RNA-seq",])

labls <- c(paste("Neither n=", nneither, sep=""),
            paste("Both n=", nboth, sep=""),
	    paste("Microarray n=", narray, sep=""),
	    paste("RNA-seq n=", nrnaseq, sep=""))

# plot
cor1 <- cor(dat$log2FoldChange, dat$maxf, use="pairwise.complete.obs")
p1 <- ggplot(dat, aes(x=log2FoldChange, y=maxf, colour=label, alpha=0.4))
p2 <- p1 + geom_point(pch=18)
p3 <- p2 + theme_bw() + scale_colour_manual(values=c("grey", "purple", "red", "blue"), labels=labls)
p4 <- p3 + xlim(c(-7,10)) + ylim(c(-7,10))
p5 <- p4 + annotate(geom="text", -5, 5, label=paste("r=", round(cor1,2), sep=""))
p6 <- p5 + geom_smooth(data=dat, aes(x=log2FoldChange, y=maxf), method="lm", inherit.aes=F)
ggsave("scatter_fold_salmonella.pdf")