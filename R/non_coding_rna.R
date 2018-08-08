##################################################
##################################################
##################################################
# Assessment of non-coding RNAs on microarray and
# RNA-seq data
##################################################
##################################################
##################################################

library(RSQLite)
library(ggplot2)

############################
# Microarrays
############################

# read ensembl to illumina ID mapping
# note these are only genes that have an associated illumina id on the array
ensembl2illumina <- read.csv("ensembl2illuminaHT12_v4.txt", header=T, stringsAsFactors=F, sep="\t")

# make a list of non-coding definitions from here
non_coding_rna <- c("non_coding", "bidirectional_promoter_lncRNA", "lincRNA", "3prime_overlapping_ncRNA", "vaultRNA")

# get proportion of each biotype
total <- nrow(ensembl2illumina)
ntable <- data.frame(table(ensembl2illumina$Gene.type))
ptable <- ntable
ptable$Freq <- round((ptable$Freq/total)*100,2)

############################
# Ensembl genes 91
############################

database <- "/gfs/mirror/annotations/hg38_ensembl91/csvdb"
sqlite <- dbDriver("SQLite")
db <- dbConnect(sqlite, database)
ens2biotype <- dbGetQuery(db, 'SELECT DISTINCT gene_id, gene_biotype FROM gene_info')
rownames(ens2biotype) <- ens2biotype$gene_id

total2 <- nrow(ens2biotype)
ntable2 <- data.frame(table(ens2biotype$gene_biotype))
ptable2 <- ntable2
ptable2$Freq <- round((ptable2$Freq/total2)*100,2)


############################
# do the plot comparison
############################

ptable$type <- "Illumina_HT12_V4"
ptable2$type <- "Ensembl_genes_91"
dat <- data.frame(rbind(ptable, ptable2))

# just look at features with > 1% frequency
dat <- dat[dat$Freq >= 1,]

p1 <- ggplot(dat, aes(x=type, y=Freq, fill=Var1))
p2 <- p1 + geom_bar(stat="identity") + theme_bw()
p2
ggsave("proportion_gene_types.pdf")


############################
# distribution of gene types
# among differentially
# expressed genes
############################

getSig <- function(diff.dat){

       sig.ids <- diff.dat$gene_id[diff.dat$padj < 0.05 & !(is.na(diff.dat$padj))]
       return(sig.ids)
       }


# cpz vs. control
cpz.control <- read.csv("cpz_vs_control.diff", header=T, stringsAsFactors=F, sep="\t")
sig.cpz <- getSig(cpz.control)
sig.cpz <- data.frame("gene_id" = sig.cpz, "gene_biotype" = ens2biotype[sig.cpz,]$gene_biotype)
sig.cpz.t <- data.frame(table(sig.cpz$gene_biotype))
sig.cpz.p <- sig.cpz.t
sig.cpz.p$Freq <- (sig.cpz.p$Freq/nrow(sig.cpz))*100


# salmonella vs. control
sal.control <- read.csv("sal_vs_control.diff", header=T, stringsAsFactors=F, sep="\t")
sig.sal <- getSig(sal.control)
sig.sal <- data.frame("gene_id" = sig.sal, "gene_biotype" = ens2biotype[sig.sal,]$gene_biotype)
sig.sal.t <- data.frame(table(sig.sal$gene_biotype))
sig.sal.p <- sig.sal.t
sig.sal.p$Freq <- (sig.sal.p$Freq/nrow(sig.sal))*100


# combine
sig.cpz.p$contrast <- "CPZ_vs_control"
sig.sal.p$contrast <- "Salmonella_vs_control"
dat.sig.p <- data.frame(rbind(sig.cpz.p, sig.sal.p))

p1 <- ggplot(dat.sig.p, aes(x=contrast, y=Freq, fill=Var1))
p2 <- p1 + geom_bar(stat="identity") + theme_bw()
p2
ggsave("proportion_gene_types_differential.pdf")


##############################################
##############################################
# fold change distributions
##############################################
##############################################

cpz.control$contrast <- "CPZ_vs_control"
sal.control$contrast <- "Salmonella_vs_control"
dat.f <- data.frame(rbind(cpz.control, sal.control))
p3 <- ggplot(dat.f, aes(x=log2FoldChange, fill=contrast, alpha=0.5))
p4 <- p3 + geom_density() + theme_bw()
p4
ggsave("fold_change_distributions.pdf")
