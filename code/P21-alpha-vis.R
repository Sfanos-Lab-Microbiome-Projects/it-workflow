rm(list=ls()) #reset of all commands in system
library(RColorBrewer) #install these packages first
library(ggplot2)
library(grid)
library(reshape2)
require(gtools)
library(ggpubr)
# ---------------------------------------------------
# configure per user ---
#setwd("it-workflow/code") # git repository already cloned to desktop
setwd(paste0(path.expand("~"),"/Desktop/it-workflow/code"))

analysisdir = "../analysis/P21-alpha-vis" #create the analysis directory you want outputs to go into
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
pdfdir = paste(analysisdir, "/pdfs/", sep="")
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)
# ---------------------------------------------------
tablefile = "../analysis/P08-filtered-sum-alpha/alpha-diversity-per-sample-region.txt"
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

  mycolors$SampleID["atcc1"] = "#FFDFD3"
  mycolors$SampleID["atcc2"] = "#FEC8D8"
  mycolors$SampleID["atcc3"] = "#E0BBE4"
  mycolors$SampleID["atcc4"] = "#D291BC"
  mycolors$SampleID["atcc5"] = "#957DAD"

  # END color scheme for major sample metadata features of interest ---------
  # --------------------------------------------------------

  outfile1 = paste(pdfdir, "alpha-diversity-by-region.1.pdf", sep="")
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

    # comparisons for statistical analysis
    my_comparisons = list(
                     c("v2", "v3" ),
                     c("v2", "v4" ),
                     c("v2", "v67"),
                     c("v2", "v8" ),
                     c("v2", "v9" ),
                     c("v3", "v4" ),
                     c("v3", "v67"),
                     c("v3", "v8" ),
                     c("v3", "v9" ),
                     c("v4", "v67"),
                     c("v4", "v8" ),
                     c("v4", "v9" ),
                     c("v67", "v8"),
                     c("v67", "v9"),
                     c("v8", "v9"))

    # begin statistical comparisons of alpha diversity estimators by region
    statRes = c()
    for (compi in 1:length(my_comparisons)){
      resultsRow = c(my_comparisons[[compi]][1], my_comparisons[[compi]][2])
      for (myFeature in unique(meltA$Measure)){
        mwP = wilcox.test(meltA[meltA$Region==my_comparisons[[compi]][1] &
                                meltA$Measure==myFeature, "Value"],
                          meltA[meltA$Region==my_comparisons[[compi]][2] &
                                meltA$Measure==myFeature, "Value"])$p.value
        ttP = t.test(meltA[meltA$Region==my_comparisons[[compi]][1] &
                           meltA$Measure==myFeature, "Value"],
                     meltA[meltA$Region==my_comparisons[[compi]][2] &
                           meltA$Measure==myFeature, "Value"])$p.value
        mn1 = mean(meltA[meltA$Region==my_comparisons[[compi]][1] &
                                meltA$Measure==myFeature, "Value"])
        mn2 = mean(meltA[meltA$Region==my_comparisons[[compi]][2] &
                                meltA$Measure==myFeature, "Value"])
        resultsRow = c(resultsRow, mn1, mn2, mwP, ttP)
      }
      statRes = rbind(statRes, resultsRow)
    }
    colnames(statRes) = c("Region 1", "Region 2", "Evenness.Mean1",      "Evenness.Mean2",     "Evenness.MWPval", "Evenness.TTPval",
                                                  "Faith's-PD.Mean1",    "Faith's-PD.Mean2",    "Faith's-PD.MWPval", "Faith's-PD.TTPval",
                                                  "Observed-OTUs.Mean1", "Observed-OTUs.Mean2", "Observed-OTUs.MWPval", "Observed-OTUs.TTPval",
                                                  "Shannon.Mean1",       "Shannon.Mean2",        "Shannon.MWPval", "Shannon.TTPval")
    write.table(statRes, file=paste(analysisdir, "/alpha-diversity-stats.csv", sep=""), col.names=TRUE, row.names=FALSE, sep=",")
