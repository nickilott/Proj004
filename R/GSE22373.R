################################################
################################################
################################################
# This analysis aims to evaluate the expression
# of neurotransmitter receptor genes in this
# publically available data set
################################################
################################################
################################################

library("biomaRt")
library("ggplot2")
library("ggrepel")
library("gcrma")

smap <- data.frame(geo=c("GSM556646", "GSM556647", "GSM556648", "GSM556663", "GSM556664", "GSM556665"),
                   cond=c("Macrophage.control", "Monocyte.control", "Macrophage.control", "Monocyte.control", "Macrophage.stim", "Monocyte.stim"),
		   rep=c("R1", "R1", "R2", "R2", "R1", "R1"))

# read matrix file
mat <- ReadAffy()

# normalise
mat <- rma(mat)
mat <- exprs(mat)
probes <- rownames(mat)
colnames(mat) <- gsub(".CEL", "", colnames(mat))

# read neurotransmitter file
ntrs <- read.csv("neurotransmitter_receptor_symbols.tsv", header=T, stringsAsFactors=F, sep="\t")

# use biomaRt to get gene symbols
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

probe2gene <- getBM(attributes = c("affy_hg_u133_plus_2", "external_gene_name"), filters = "affy_hg_u133_plus_2", values = probes, mart = mart)
colnames(probe2gene) <- c("probe", "gene")
probe2gene$gene <- toupper(probe2gene$gene)
write.table(probe2gene, file = "probe2gene.tsv", sep = "\t", row.names = F)

# plot the distribution of (mean) expression values
mat <- data.frame(mat)
mexprs <- data.frame("average" = rowMeans(mat))
mexprs$type <- "Other"
mexprs$probe <- rownames(mat)

# get probes for NTRs
ntr.probes <- probe2gene$probe[probe2gene$gene %in% ntrs$NTR]
mexprs$type <- ifelse(mexprs$probe %in% ntr.probes, "NTR", mexprs$type)

# annotate NTR gene names where expression > 4
probe2gene.ntr <- probe2gene[probe2gene$probe %in% ntr.probes,]

# some weird inclusion of HTR7P1 - remove
probe2gene.ntr <- probe2gene.ntr[grep("HTR7P1", probe2gene.ntr$gene, invert=T),]

rownames(probe2gene.ntr) <- probe2gene.ntr$probe
mexprs$gene <- ifelse(mexprs$type=="NTR" & mexprs$average >=4, probe2gene.ntr[mexprs$probe,]$gene , NA)

p1 <- ggplot(mexprs, aes(x=average, fill=type, alpha=0.5))
p2 <- p1 + geom_density() + theme_bw() + scale_fill_manual(values=c("red", "grey"))
p3 <- p2 + geom_text_repel(aes(x=average, y=0, label=gene))
p4 <- p3 + xlab("Average RMA normalised log2(expression)")
p5 <- p4 + ylab("Density")
p6 <- p5 + ggtitle("Expression of neurotransmitter receptors in human gut macrophages/monocytes (GSE22373)")
p6
ggsave("GSE22373_colonic_macrophages.pdf")


#########################################
#########################################
#########################################
# compare monocytes to macrophages
#########################################
#########################################
#########################################

mono <- mat[,as.character(smap$geo[grep("Monocyte", smap$cond)])]
mac <- mat[,as.character(smap$geo[grep("Macrophage", smap$cond)])]

mono.ave <- data.frame("Cell.type" = "Monocyte", "Monocyte" = rowMeans(mono))
mac.ave <- data.frame("Cell.type" = "Macrophage", "Macrophage" = rowMeans(mac))

dat <- data.frame(cbind(mono.ave, mac.ave))

dat$type <- ifelse(rownames(dat) %in% probe2gene.ntr$probe, "NTR", "Other")
dat$gene <- ifelse(dat$type == "NTR" & dat$Monocyte >= 4 | dat$Macrophage >= 4, probe2gene.ntr[rownames(dat),]$gene, NA)

p1 <- ggplot(dat, aes(x=Monocyte, y=Macrophage, colour=type, size=type, alpha=type))
p2 <- p1 + geom_point(pch=18)
p3 <- p2 + scale_size_manual(values=c(3,1))
p4 <- p3 + scale_colour_manual(values=c("red4", "grey"))
p5 <- p4 + theme_bw()
p6 <- p5 + scale_alpha_manual(values=c(0.5, 0.1))
p7 <- p6 + geom_text_repel(aes(label=gene))
p8 <- p7 + xlab("RMA normalised expression log2 (Monocyte)")
p9 <- p8 + ylab("RMA normalised expression log2 (Macrophage)")
ggsave("GSE22373_monocyte_vs_macrophage.pdf")