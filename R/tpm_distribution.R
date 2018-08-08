################################################
################################################
################################################
# This analysis aims to evaluate the expression
# of neurotransmitter receptor genes in MDMs
# RNA-seq data
################################################
################################################
################################################

library("ggplot2")
library("ggrepel")
library("RSQLite")

# read matrix file
mat <- read.csv("genes.tsv.gz", header=T, stringsAsFactors=F, sep="\t")

# get genes associated with ensembl gene ids
database <- "/gfs/mirror/annotations/hg38_ensembl91/csvdb"
sqlite <- dbDriver("SQLite")
db <- dbConnect(sqlite, database)
ens2name <- dbGetQuery(db, 'SELECT DISTINCT gene_id, gene_name FROM gene_info')

# read neurotransmitter file
ntrs <- read.csv("neurotransmitter_receptor_symbols.tsv", header=T, stringsAsFactors=F, sep="\t")

rownames(mat) <- mat$id
mat <- mat[,2:ncol(mat)]
mat <- log10(mat+1)

mexprs <- data.frame("average" = rowMeans(mat))
mexprs$type <- "Other"
mexprs$gene_id <- rownames(mat)

# get ensembl for NTRs
ntr.ens <- ens2name$gene_id[ens2name$gene_name %in% ntrs$NTR]
mexprs$type <- ifelse(mexprs$gene_id %in% ntr.ens, "NTR", mexprs$type)

# annotate NTR gene names for NTRs > log10(TPM) >=0.5
ens2name.ntr <- ens2name[ens2name$gene_id %in% ntr.ens,]
rownames(ens2name.ntr) <- ens2name.ntr$gene_id
mexprs$gene <- ifelse(mexprs$type=="NTR" & mexprs$average >= 0.5, ens2name.ntr[mexprs$gene_id,]$gene_name , NA)

# add density for annotation purposes
dx <- density(mexprs$average)
ds <- approx(dx$x,dx$y, xout=mexprs$average)
mexprs$ds <- ds$y
mexprs$ds <- ifelse(mexprs$type=="NTR" & mexprs$average >= 0.5, mexprs$ds, NA)

p1 <- ggplot(mexprs, aes(x=average, fill=type, alpha=0.5))
p2 <- p1 + geom_density() + theme_bw() + scale_fill_manual(values=c("red", "grey"))
p3 <- p2 + geom_text_repel(aes(x=average, y=0, label=gene))
p4 <- p3 + xlab("Average transcripts per million (log10(TPM + 1)) ")
p5 <- p4 + ylab("Density")
p6 <- p5 + ggtitle("Expression of neurotransmitter receptors in MDMs")
p6
ggsave("mdm_tpms.pdf")
