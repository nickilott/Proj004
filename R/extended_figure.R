##########################################################
##########################################################
##########################################################
# A subset of figures used for the Extended figure for
# Sumeet's revision
##########################################################
##########################################################
##########################################################

library(ggplot2)
library(gplots)
library(VennDiagram)
source("/gfs/devel/nilott/NGSKit/R/deseq2_helper.R")

# read data
mat <- read.csv("genes_rlog.tsv", header=T, stringsAsFactors=F, sep="\t", row.names=17)

##################################
##################################
# principle components analysis
##################################
##################################

pc <- PCA(mat)

# get donor as covariate
donor <- unlist(strsplit(rownames(pc$x), "\\."))
donor <- donor[seq(3, length(donor), 3)]

# plot PCA
pl1 <- plotPCAWithCovariate(pc, donor, pcs=c("PC1", "PC2"))
pl2 <- pl1 + scale_shape_manual(values=c(3,15, 17,18))
pl3 <- pl2 + scale_colour_manual(values=c("blue", "blue4", "red", "red4"))
ggsave("PC1_PC2.pdf")

pl4 <- plotPCAWithCovariate(pc, donor, pcs=c("PC1", "PC3"))
pl5 <- pl4 + scale_shape_manual(values=c(3,15, 17,18))
pl6 <- pl5 + scale_colour_manual(values=c("blue", "blue4", "red", "red4"))
ggsave("PC1_PC3.pdf")

pl7 <- plotPCAWithCovariate(pc, donor, pcs=c("PC2", "PC3"))
pl8 <- pl7 + scale_shape_manual(values=c(3,15, 17,18))
pl9 <- pl8 + scale_colour_manual(values=c("blue", "blue4", "red", "red4"))
ggsave("PC2_PC3.pdf")


##################################
##################################
# Heatmap
##################################
##################################

readData <- function(x){

	 dat <- read.csv(x, header=T, stringsAsFactors=F, sep="\t")
	 return(dat)
	 }

# read annotated matrix
mat.annotated <- readData("genes_rlog_annotated.tsv")

# read diff sets
cpz.control <- readData("cpz_vs_control.diff")
sal.control <- readData("sal_vs_control.diff")
salcpz.sal <- readData("salcpz_vs_sal.diff")
salcpz.cpz <- readData("salcpz_vs_cpz.diff")

# get union of differentially regulated genes
un <- union(union(union(cpz.control$gene_name, sal.control$gene_name), salcpz.sal), salcpz.cpz)

# subset matrix
mat.annotated <- mat.annotated[which(mat.annotated$gene_name %in% un),]

# reform
rownames(mat.annotated) <- mat.annotated$gene_id
mat.annotated <- mat.annotated[,1:ncol(mat.annotated)-1]
mat.annotated <- mat.annotated[,1:ncol(mat.annotated)-1]

# plot
col.anno <- makeConds(colnames(mat.annotated))
col.anno <- ifelse(col.anno == "Control", "blue", col.anno)
col.anno <- ifelse(col.anno == "cpz", "red", col.anno)
col.anno <- ifelse(col.anno == "Salmonella", "blue4", col.anno)
col.anno <- ifelse(col.anno == "cpzSalmonella", "red4", col.anno)

# colours for values
colours <- colorRampPalette(c("blue", "white", "red"))(75)

# scale data - zscore
toplot <- data.frame(t(apply(mat.annotated, 1, scale)))
colnames(toplot) <- colnames(mat.annotated)

# heatmap
png("heatmap.png", height=12, width=10, units="in", res=300)
heatmap.2(as.matrix(toplot), trace="none", col=colours, margins=c(15,15), scale="none", ColSideColors=col.anno, labRow=F)

# add legend
legend("topright", legend=c("Control", "CPZ", "Salmonella", "Salmonella + CPZ"), fill=c("blue", "red", "blue4", "red4"))
dev.off()

################################
################################
# venn diagrams
################################
################################

x <- list("CPZ vs Control" = cpz.control$gene_name, "Salmonella vs. Control" = sal.control$gene_name)
y <- list("Salmonella + CPZ vs CPZ" = salcpz.cpz$gene_name, "Salmonella + CPZ vs. Salmonella" = salcpz.sal$gene_name)

venn.diagram(x, filename="venn_vs_control.tiff")
venn.diagram(y, filename="venn_vs_salmonella.tiff")

