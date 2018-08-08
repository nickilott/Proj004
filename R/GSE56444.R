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

# read matrix file
mat <- read.csv("sample_probe_profile.matrix", header=T, stringsAsFactors=F, sep="\t")

probes <- mat$ID_REF

# read neurotransmitter file
ntrs <- read.csv("neurotransmitter_receptor_symbols.tsv", header=T, stringsAsFactors=F, sep="\t")


# use biomaRt to get gene symbols
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

probe2gene <- getBM(attributes = c("affy_mogene_1_0_st_v1", "external_gene_name"), filters = "affy_mogene_1_0_st_v1", values = probes, mart = mart)
colnames(probe2gene) <- c("probe", "gene")
probe2gene$gene <- toupper(probe2gene$gene)
write.table(probe2gene, file = "probe2gene.tsv", sep = "\t", row.names = F)

# plot the distribution of (mean) expression values
rownames(mat) <- mat$ID_REF
mat <- mat[,2:ncol(mat)]
mexprs <- data.frame("average" = rowMeans(mat))
mexprs$type <- "Other"
mexprs$probe <- rownames(mat)

# get probes for NTRs
ntr.probes <- probe2gene$probe[probe2gene$gene %in% ntrs$NTR]
mexprs$type <- ifelse(mexprs$probe %in% ntr.probes, "NTR", mexprs$type)

# annotate NTR gene names where expression > 5
probe2gene.ntr <- probe2gene[probe2gene$probe %in% ntr.probes,]
rownames(probe2gene.ntr) <- probe2gene.ntr$probe
mexprs$gene <- ifelse(mexprs$type=="NTR" & mexprs$average >=5, probe2gene.ntr[mexprs$probe,]$gene , NA)

# add density for annotation purposes
dx <- density(mexprs$average)
ds <- approx(dx$x,dx$y, xout=mexprs$average)
mexprs$ds <- ds$y
mexprs$ds <- ifelse(mexprs$type=="NTR" & mexprs$average >=5, mexprs$ds, NA)


p1 <- ggplot(mexprs, aes(x=average, fill=type, alpha=0.5))
p2 <- p1 + geom_density() + theme_bw() + scale_fill_manual(values=c("red", "grey"))
p3 <- p2 + geom_text_repel(aes(x=average, y=0, label=gene))
p4 <- p3 + xlab("Average RMA normalised log2(expression)")
p5 <- p4 + ylab("Density")
p6 <- p5 + ggtitle("Expression of neurotransmitter receptors in sorted colonic murine macrophages")
p6
ggsave("GSE56444_colonic_macrophages.pdf")