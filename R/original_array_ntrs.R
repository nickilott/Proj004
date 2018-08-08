####################################################
####################################################
# Analysis of the expression of neurotransmitter
# receptors in the original array data
####################################################
####################################################

library("biomaRt")
library("ggplot2")
library("ggrepel")

# read matrix file
mat <- read.csv("sample_probe_profile.matrix", header=T, stringsAsFactors=F, sep="\t")
rownames(mat) <- mat$ids
mat <- mat[,1:ncol(mat)-1]

# read neurotransmitter file
ntrs <- read.csv("neurotransmitter_receptor_symbols.tsv", header=T, stringsAsFactors=F, sep="\t")

# read probe2gene annotation
probe2gene <- read.csv("probe2gene.map", header=T, stringsAsFactors=F, sep="\t")
rownames(probe2gene) <- probe2gene$probe

# get neurotransmitter receptor probes
ntr.probes <- probe2gene$probe[probe2gene$gene %in% ntrs$NTR]

# build average expression and annotate NTRs
mexprs <- data.frame("average" = rowMeans(mat))
mexprs$type <- "Other"
mexprs$probe <- rownames(mat)
mexprs$type <- ifelse(mexprs$probe %in% ntr.probes, "NTR", mexprs$type)

# annotate NTR gene names where expression > 7
probe2gene.ntr <- probe2gene[probe2gene$probe %in% ntr.probes,]
rownames(probe2gene.ntr) <- probe2gene.ntr$probe
mexprs$gene <- ifelse(mexprs$type=="NTR" & mexprs$average >=7, probe2gene.ntr[mexprs$probe,]$gene , NA)

# plot
p1 <- ggplot(mexprs, aes(x=average, fill=type, alpha=0.5))
p2 <- p1 + geom_density() + theme_bw() + scale_fill_manual(values=c("red", "grey"))
p3 <- p2 + geom_text_repel(aes(x=average, y=0, label=gene))
p4 <- p3 + xlab("Average RMA normalised log2(expression)")
p5 <- p4 + ylab("Density")
p6 <- p5 + ggtitle("Expression of neurotransmitter receptors in MDMs (Illumina HT12 V4 array)")
p6
ggsave("Illumina_HT12_v4_array.pdf")
