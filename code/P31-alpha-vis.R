rm(list=ls()) #reset of all commands in system
library(RColorBrewer) #install these packages first
library(ggplot2)
library(grid)
library(reshape2)
# ---------------------------------------------------
# configure per user ---
#setwd("it-workflow/code") # git repository already cloned to desktop
analysisdir = "../analysis/P31-alpha-vis" #create the analysis directory you want outputs to go into
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
pdfdir = paste(analysisdir, "/pdfs/", sep="")
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)
# ---------------------------------------------------
tablefile = "../analysis/P18-sum-alpha/alpha-diversity-per-sample-region.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
meltA <- melt(A, id.vars=c("Region", "SampleID"))
colnames(meltA) = c("Region", "SampleID", "Measure", "Value")
meltA$Measure = gsub("exported-","", meltA$Measure)
meltA$Measure = gsub("-vector","", meltA$Measure)
meltA$Measure = gsub("evenness","Evenness", meltA$Measure)
meltA$Measure = gsub("faith-pd","Faith's-PD", meltA$Measure)
meltA$Measure = gsub("observed-otus","Observed-OTUs", meltA$Measure)
meltA$Measure = gsub("shannon","Shannon", meltA$Measure)
  # --------------------------------------------------------
  # BEGIN color scheme for major sample metadata features of interest ---------
  mycolors = c()
  mycolors$Region["v2"]  = "#FF9AA2"
  mycolors$Region["v3"]  = "#FFDAC1"
  mycolors$Region["v4"]  = "#FFF7C1"
  mycolors$Region["v67"] = "#E2F0CB"
  mycolors$Region["v8"]  = "#B5EAD7"
  mycolors$Region["v9"]  = "#C7CEEA"

  mycolors$SampleID["Patient1-fresh"] = "#FFA762"
  mycolors$SampleID["Patient1-frozen"] = "#FFD2AF"
  mycolors$SampleID["Patient2-fresh"] = "#6D5FA6"
  mycolors$SampleID["Patient2-frozen"] = "#B3ACD0"
  mycolors$SampleID["Patient3-fresh"] = "#518FBF"
  mycolors$SampleID["Patient3-frozen"] = "#62FFF5"

  # END color scheme for major sample metadata features of interest ---------
  # --------------------------------------------------------
  outfile1 = paste(pdfdir, "alpha-diversity-by-region.pdf", sep="")
  p1 <- ggplot(meltA, aes(x=Region, y=Value)) +
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
    ylab("Measure Value") +
    facet_wrap(~Measure, ncol=2, scales="free_y")
  ggsave(outfile1, plot=p1)
