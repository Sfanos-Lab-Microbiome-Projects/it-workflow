rm(list=ls()) #reset of all commands in system
library(RColorBrewer) #install these packages first
library(ggplot2)
library(grid)
library(reshape2)
# ---------------------------------------------------
# configure per user ---
#setwd("it-workflow/code") # git repository already cloned to desktop
analysisdir = "../analysis/P23-filtered-tax-class" #create the analysis directory you want outputs to go into
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
pdfdir = paste(analysisdir, "/pdfs/", sep="")
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)
# ---------------------------------------------------
tablefile = "../analysis/P07-filtered-sum-taxa-db4.0/filtered-species-for-r.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
meltA <- melt(A, id.vars=c("Region", "SampleID"))
colnames(meltA) = c("Region", "SampleID", "Taxa", "PercentAbundance")
# --------------------------------------------------------
# load meta-data ---
metadatafile                = "../analysis/P07-filtered-sum-taxa-db4.0/filtered-taxa-meta.txt"
meta                        = read.table(metadatafile, header=TRUE, sep="\t")
meta                        = data.frame(meta)
# ---------------------------------------------------
# match the meta data frame and meltA to get the needed variables for plotting
mergeAwmeta = merge(meltA,meta,by="Taxa")
# order the factor levels for visualization --
mergeAwmeta$Taxa      = factor(mergeAwmeta$Taxa)
mergeAwmeta$Result     = factor(mergeAwmeta$Result,   levels=c("TP","FN","FP"))
mergeAwmeta$Region     = factor(mergeAwmeta$Region,   levels=c("v2","v3","v4","v67","v8","v9","Expected"))
mergeAwmeta$SampleID    = factor(mergeAwmeta$SampleID,  levels=c("atcc1","atcc2","atcc3","atcc4","atcc5","Expected"))

Awmeta = mergeAwmeta[mergeAwmeta$Result != "FP" & mergeAwmeta$Result != "FN", ]
# BEGIN color scheme for major sample metadata features of interest ---------
mycolors = c()
mycolors$Region["Expected"]  = "#FF9AA2"
mycolors$Region["v2"]  = "#FF9AA2"
mycolors$Region["v3"]  = "#FFDAC1"
mycolors$Region["v4"]  = "#FFF7C1"
mycolors$Region["v67"] = "#E2F0CB"
mycolors$Region["v8"]  = "#B5EAD7"
mycolors$Region["v9"]  = "#C7CEEA"

mycolors$SampleID["atcc1"] = "#F2978B"
mycolors$SampleID["atcc2"] = "#D954A0"
mycolors$SampleID["atcc3"] = "#6D5FA6"
mycolors$SampleID["atcc4"] = "#518FBF"
mycolors$SampleID["atcc5"] = "#63BBBF"

# END color scheme for major sample metadata features of interest ---------

# --------------------------------------------------------
outfile1 = paste(pdfdir, "percent-taxa-per-region.pdf", sep="")
p1 <- ggplot(Awmeta, aes(x=Region, y=PercentAbundance)) +
  geom_boxplot(mapping=aes(color=Region), alpha=1, outlier.size = NA, coef=1000) +
  geom_point(aes(color=Region), alpha=1, size=2) +
  theme_bw() +
  scale_color_manual(values=mycolors$Region) +
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y  = element_text(size=11, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        plot.title   = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1.0) +
  xlab(NULL) +
  ylab(NULL) 
ggsave(outfile1, plot=p1)

