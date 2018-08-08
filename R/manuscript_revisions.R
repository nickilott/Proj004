#################################################
#################################################
#################################################
# Some questions from reviewers were raised
# and these are addressed using this script
#################################################
#################################################
#################################################


readData <- function(infile){

	 dat <- read.csv(infile, header=T, stringsAsFactors=F, sep="\t")
	 return(dat)
	 }


neurotransmitters <- c("HTR7",
		       "HTR6",
		       "HTR5A",
		       "HTR4",
		       "HTR3C",
		       "HTR3B",
		       "HTR3A",
		       "HTR2C",
		       "HTR2B",
		       "HTR2A",
		       "HTR1F",
		       "HTR1E",
		       "HTR1D",
		       "HTR1C",
		       "HTR1B",
		       "HTR1A",
		       "HRH4",
		       "HRH3",
		       "HRH2",
		       "HRH1",
		       "DRD5",
		       "DRD4",
		       "DRD3",
		       "DRD2",
		       "DRD1",
		       "ADRB1",
		       "ADRA2C",
		       "ADRA2B",
		       "ADRA2A",
		       "ADRA1D",
		       "ADRA1B",
		       "ADRA1A")

#################################################
# MA plot of Samonella induced changes in
# untreated macs
#################################################

results.file = "/gfs/work/nilott/proj004/analysis/differential_expression1/differential_expression.dir/UntreatedInfected-UntreatedUnInfected.result"

dat <- readData(results.file)
dat$significant <- ifelse(dat$adj.P.Val < 0.05, "Yes", "No")

nup <- nrow(dat[dat$significant=="Yes" & dat$logFC > 0,])
ndown <- nrow(dat[dat$significant=="Yes" & dat$logFC < 0,])

nup <- paste("Upregulated = ", nup, sep="")
ndown <- paste("Downregulated = ", ndown, sep="")

library(ggplot2)

plot1 <- ggplot(dat, aes(x=AveExpr, y=logFC, colour=significant))
plot2 <- plot1 + geom_point(pch=18, size=2)
plot3 <- plot2 + scale_colour_manual(values=c("grey", "darkRed"))
plot4 <- plot3 + theme_bw()
plot5 <- plot4 + theme(text=element_text(size=10))
plot6 <- plot5 + geom_hline(yintercept=c(-1,1,0), linetype="dashed")
plot7 <- plot6 + ylim(c(-7,7)) + xlab("Mean expression") + ylab("Log2 fold change in Salmonella infected vs. Uninfected")
plot8 <- plot7 + annotate("text", x=9, y=6, label=nup, size=3)
plot9 <- plot8 + annotate("text", x=9, y=-6, label=ndown, size=3)

# write out differentially expressed genes
diff <- dat[dat$significant == "Yes",]
diff <- diff[order(diff$adj.P.Val, decreasing=F),]
write.table(diff, file="salmonella_regulated.tsv", quote=F, sep="\t", row.names=F)


#################################################
# MA plot of Samonella induced changes in
# untreated of neurotransmitters and their receptors
#################################################

library(ggrepel)

dat.n <- dat[dat$gene %in% neurotransmitters,]
p1 <- ggplot(dat.n, aes(x=AveExpr, y=logFC, colour=significant))
p2 <- p1 + geom_point(pch=18, size=2)
p3 <- p2 + scale_colour_manual(values=c("grey", "darkRed"))
p4 <- p3 + theme_bw()
p5 <- p4 + theme(text=element_text(size=10))
p6 <- p5 + geom_hline(yintercept=c(-1,1,0), linetype="dashed")
p7 <- p6 + xlab("Mean expression") + ylab("Log2 fold change in Salmonella infected vs. Uninfected")
p8 <- p7 + geom_text_repel(aes(label=gene))
ggsave("neurotransmitters.pdf")


#################################################
# Pathways analysis of Salmonella infected
# untreated macs
#################################################

pathway.results <- "/gfs/work/nilott/proj004/analysis/differential_expression1/pathways.dir/UntreatedInfected_UntreatedUnInfected_result_up.h.all.v5.1.symbols.results"
p.dat <- readData(pathway.results)
p.dat <- p.dat[p.dat$code=="+",]
p.dat <- p.dat[order(p.dat$ratio, decreasing=T),]
p.dat$label <- paste(round(p.dat$fdr, 3), " (", sep="")
p.dat$label <- paste(p.dat$label, p.dat$scount, sep="")
p.dat$label <- paste(p.dat$label, ")", sep="")

p.dat$goid <- factor(p.dat$goid, levels=p.dat$goid)

plot10 <- ggplot(p.dat, aes(x=goid, y=ratio, label=scount))
plot11 <- plot10 + geom_bar(stat="identity") + geom_text(vjust=-1)
plot12 <- plot11 + theme_bw() + theme(text=element_text(angle=90))
plot13 <- plot12 + ylim(c(0,5))
plot14 <- plot13 + xlab("") + ylab("Fold enrichment")

##################################################
# Compare responses between CPZ and untreated
# when infected with Salmonella
##################################################

cpz <- "/gfs/work/nilott/proj004/analysis/differential_expression1/differential_expression.dir/CPZUninfected-CPZInfected.result"
untreated <- "/gfs/work/nilott/proj004/analysis/differential_expression1/differential_expression.dir/UntreatedInfected-UntreatedUnInfected.result"

cpz <- readData(cpz)
untreated <- readData(untreated)

# need to flip fold changes in CPZ results file
cpz$logFC <- cpz$logFC*-1

# combine data
dat <- data.frame("untreated" = untreated$logFC,
                 "CPZ" = cpz$logFC)

r.squared <- round(cor(dat$untreated, dat$CPZ)^2,3)
r.squared <- paste("r2 = ", r.squared, sep="" )

plot15 <- ggplot(dat, aes(x=untreated, y=CPZ))
plot16 <- plot15 + geom_point(pch=18)
plot17 <- plot16 + theme_bw()
plot18 <- plot17 + theme(text=element_text(size=10))
plot19 <- plot18 + annotate("text", x=-1, y=6, label=r.squared)
plot20 <- plot19 + stat_smooth(method="lm")
plot21 <- plot20 + xlab("Log2 Fold change in Salmonella infected vs. Uninfected")
plot22 <- plot21 + ylab("Log2 Fold change in Salmonella infected vs. Uninfected (CPZ treatment)")


#################################################
# put the plots together
#################################################

library(gridExtra)

pdf("revision_plots.pdf", height=5, width=15)
grid.arrange(plot9, plot14, plot22, ncol=3, nrow=1)
dev.off()


################################################
# Expression of neurotransmitter receptors
################################################








