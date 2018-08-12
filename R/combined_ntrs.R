##########################################################
##########################################################
##########################################################
# analysis to produce plots for the expression level
# analysis of neurotransmitter receptors
##########################################################
##########################################################
##########################################################

library("biomaRt")
library("ggplot2")
library("ggrepel")
library("gplots")

# common to all analyses
# read neurotransmitter file
ntrs <- read.csv("neurotransmitter_receptor_symbols.tsv", header=T, stringsAsFactors=F, sep="\t")

########################################
# function to get gene names from
# matrix with probe ids as rownames

getGenes <- function(probes, probe2gene){
	    genes <- c()
            for (i in probes){
	      gene <- probe2gene$gene[probe2gene$probe == i]
	      genes <- append(genes, gene)
				        }
	     return(genes)
	     }
							      
##########################################
##########################################
# murine macrophages - GSE56444
##########################################
##########################################

mat <- read.csv("murine_gut_macs.matrix", header=T, stringsAsFactors=F, sep="\t", row.names=1)

probe2gene <- read.csv("GSE56444_probe2gene.map", header=T, stringsAsFactors=F, sep="\t")

# get probes for NTRs
ntr.probes <- as.character(probe2gene$probe[probe2gene$gene %in% ntrs$NTR])

mat.ntr <- mat[ntr.probes,]

mat.ntr <- mat.ntr[order(rowSums(mat.ntr), decreasing=T),]

genes <- getGenes(rownames(mat.ntr), probe2gene)

colours <- colorRampPalette(c("blue", "white", "red"))(75)

pdf("GSE56444_heatmap_ntrs.pdf", height=15, width=10)
heatmap.2(as.matrix(mat.ntr), trace="none", scale="none", margins=c(15,15), col=colours, Rowv=F, labRow=genes)
dev.off()


##########################################
##########################################
# human macrophages/monocytes - GSE22373
##########################################
##########################################

library("gcrma")

# read matrix file
mat <- ReadAffy()

smap <- data.frame(geo=c("GSM556646", "GSM556647", "GSM556648", "GSM556663", "GSM556664", "GSM556665"),
                   cond=c("Macrophage.control", "Monocyte.control", "Macrophage.control", "Monocyte.control", "Macrophage.stim", "Monocyte.stim"),
		   rep=c("R1", "R1", "R2", "R2", "R1", "R1"))
rownames(smap) <- smap$geo

# normalise
mat <- rma(mat)
mat <- exprs(mat)
probes <- rownames(mat)
colnames(mat) <- gsub(".CEL", "", colnames(mat))

probe2gene <- read.csv("GSE22373_probe2gene.map", header=T, stringsAsFactors=F, sep="\t")

# get probes for NTRs
ntr.probes <- as.character(probe2gene$probe[probe2gene$gene %in% ntrs$NTR])

mat.ntr <- mat[ntr.probes,]

mat.ntr <- mat.ntr[order(rowSums(mat.ntr), decreasing=T),]

genes <- getGenes(rownames(mat.ntr), probe2gene)

colours <- colorRampPalette(c("blue", "white", "red"))(75)

pdf("GSE22373_heatmap_ntrs.pdf", height=15, width=10)
heatmap.2(as.matrix(mat.ntr), trace="none", scale="none", margins=c(15,15), col=colours, Rowv=F, labRow=genes, labCol=smap[colnames(mat.ntr),]$cond)
dev.off()

##########################################
##########################################
# MDM original array
##########################################
##########################################

mat <- read.csv("original_array.matrix", header=T, stringsAsFactors=F, sep="\t")

probe2gene <- read.csv("original_array_probe2gene.map", header=T, stringsAsFactors=F, sep="\t")

# get probes for NTRs
ntr.probes <- as.character(probe2gene$probe[probe2gene$gene %in% ntrs$NTR])

mat.ntr <- mat[mat$ids %in% ntr.probes,]

rownames(mat.ntr) <- mat.ntr$ids

mat.ntr <- mat.ntr[,1:ncol(mat.ntr)-1]

mat.ntr <- mat.ntr[order(rowSums(mat.ntr), decreasing=T),]

genes <- getGenes(rownames(mat.ntr), probe2gene)

colours <- colorRampPalette(c("blue", "white", "red"))(75)

pdf("original_array_mdm_heatmap_ntrs.pdf", height=15, width=10)
heatmap.2(as.matrix(mat.ntr), trace="none", scale="none", margins=c(15,15), col=colours, Rowv=F, labRow=genes)
dev.off()

##########################################
##########################################
# RNA-seq tpm
##########################################
##########################################

library(RSQLite)

mat <- read.csv("rnaseq_tpms.tsv.gz", header=T, stringsAsFactors=F, sep="\t", row.names=1)

# get genes associated with ensembl gene ids
database <- "/gfs/mirror/annotations/hg38_ensembl91/csvdb"

sqlite <- dbDriver("SQLite")

db <- dbConnect(sqlite, database)

ens2name <- dbGetQuery(db, 'SELECT DISTINCT gene_id, gene_name FROM gene_info')

# get ensembl gene_ids for NTRs
ntr.ens <- as.character(ens2name$gene_id[ens2name$gene_name %in% ntrs$NTR])

mat.ntr <- mat[rownames(mat) %in% ntr.ens,]

mat.ntr <- mat.ntr[order(rowSums(mat.ntr), decreasing=T),]

mat.ntr <- log10(mat.ntr + 1)

genes <- c()
for (i in rownames(mat.ntr)){
    gene <- ens2name$gene_name[ens2name$gene_id == i]
    genes <- append(genes, gene)
    }
    
pdf("rnaseq_mdm_heatmap_ntrs.pdf", height=15, width=10)
heatmap.2(as.matrix(mat.ntr), trace="none", scale="none", margins=c(15,15), col=colours, Rowv=F, labRow=genes)
dev.off()



