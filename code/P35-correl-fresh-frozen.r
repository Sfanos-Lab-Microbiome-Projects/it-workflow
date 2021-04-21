rm(list=ls()) #reset of all commands in system
library(RColorBrewer) #install these packages first
library(ggplot2)
library(grid)
library(reshape2)
require(GGally)
require(gtools)
library(ggpubr)
# ---------------------------------------------------
# configure per user ---
#setwd("it-workflow/code") # git repository already cloned to desktop
setwd(paste0(path.expand("~"),"/Desktop/it-workflow/code"))

analysisdir = "../analysis/P35-correl-fresh-frozen" # create the analysis directory you want outputs to go into
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
pdfdir = paste(analysisdir, "/pdfs/", sep="") #subdirectory in analysisdir
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)

txlevels   = c("alpha-diversity",
               "phylum",
               "class",
               "order",
               "family",
               "genus",
               "species",
               "otu-table")

for (txlevel in txlevels){
  # ---------------------------------------------------
  tablefile = paste("../analysis/P17-sum-taxa-db4.0.pl/final/", txlevel, ".txt", sep="")
  if (grepl("alpha", txlevel)){
    tablefile = "../analysis/P18-sum-alpha/alpha-diversity-per-sample-region.txt"
  }
  A             = read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
  Awmeta        = A
  # skip v9
  Awmeta        = Awmeta[Awmeta$Region != "v9" & Awmeta$Region != "v8", ]
  Awmeta$Region = factor(Awmeta$Region)
  Awmeta$Type   = "fresh"
  Awmeta[grepl("frozen", Awmeta$SampleID),"Type"] = "frozen"
  Awmeta$PID    = gsub("-(fresh|frozen)", "", as.character(Awmeta$SampleID))
  Awmeta$Type   = factor(Awmeta$Type)
  Awmeta        = Awmeta[order(Awmeta$Region, Awmeta$Type, Awmeta$PID),]

  meltA           = melt(Awmeta, id.vars = c("Region", "SampleID", "Type", "PID"))
  colnames(meltA) = c("Region", "SampleID", "Type", "PID", "Taxa", "Abundance")
  meltA$Taxa      = as.character(meltA$Taxa)
  meltA$Abundance = meltA$Abundance+0.001
  pairedmeltA = c()
  for (region in unique(meltA$Region)){
    for (pid in unique(meltA$PID)){
      myslice = cbind(meltA[meltA$Region == region &
                      meltA$PID  == pid &
                      meltA$Type == "fresh",],
                      meltA[meltA$Region == region &
                      meltA$PID  == pid &
                      meltA$Type == "frozen","Abundance"])
      pairedmeltA = rbind(pairedmeltA, myslice)
    }
  }
  colnames(pairedmeltA) = c("Region", "SampleID", "Type", "PID", "Taxa", "Abundance", "AbundanceFrz")
  write.table(pairedmeltA, paste(pdfdir, txlevel, ".01.txt", sep=""), sep="\t", col.names=TRUE, row.names=FALSE)



  # --------------------------------------------------------
  # BEGIN color scheme for major sample metadata features of interest ---------
  mycolors = c()
  mycolors$Type["fresh"]  = "#74c7b8"
  mycolors$Type["frozen"] = "#749fc7"

  mycolors$Region["v2"]  = "#FF9AA2"
  mycolors$Region["v3"]  = "#FFDAC1"
  mycolors$Region["v4"]  = "#FFF7C1"
  mycolors$Region["v67"] = "#E2F0CB"
  mycolors$Region["v8"]  = "#B5EAD7"
  mycolors$Region["v9"]  = "#C7CEEA"

  mycolors$PID["Patient1"] = "#98ddca"
  mycolors$PID["Patient2"] = "#ffd3b4"
  mycolors$PID["Patient3"] = "#ffaaa7"
  # END color scheme for major sample metadata features of interest ---------
  # --------------------------------------------------------
  # BEGIN VIS CODE -------------------------------------------
  # ggplot2 ---
  p1 <- ggplot(pairedmeltA[sample(1:nrow(pairedmeltA)),], aes_string(x="Abundance", y="AbundanceFrz")) +
  geom_point(aes(color=PID), alpha=0.75, size=2.5) +
  theme_bw() +
  scale_color_manual(values=mycolors$PID) +
  theme(axis.text.x  = element_text(size=8, colour="black", angle=45, hjust=1, vjust=1),
        axis.text.y  = element_text(size=8, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        plot.title   = element_text(size=9, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=1.0) +
  xlab("Fresh % Ab.") +
  ylab("Frozen % Ab.") +
  facet_wrap(~Region, ncol=2) +
  geom_abline(slope=1, intercept=0, col="black", linetype = 2) +
  scale_x_log10(breaks=c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100)) +
  scale_y_log10(breaks=c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100))

  if (grepl("alpha",txlevel)){
    pairedmeltA = pairedmeltA[grepl("observed.otus", pairedmeltA$Taxa), ]
    p1 <- ggplot(pairedmeltA, aes_string(x="Abundance", y="AbundanceFrz")) +
    geom_point(aes(color=PID), alpha=0.95, size=3) +
    theme_bw() +
    scale_color_manual(values=mycolors$PID) +
    theme(axis.text.x  = element_text(size=8, colour="black", angle=45, hjust=1, vjust=1),
          axis.text.y  = element_text(size=8, colour="black"),
          axis.title.x = element_text(size=12, colour="black"),
          axis.title.y = element_text(size=12, colour="black"),
          plot.title   = element_text(size=9, colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio=1.0) +
    xlab("Fresh % Ab.") +
    ylab("Frozen % Ab.") +
    facet_wrap(~Region, ncol=2) +
    geom_abline(slope=1, intercept=0, col="black", linetype = 2)
  }
  ggsave(paste(pdfdir, txlevel, ".01.pdf", sep=""), plot=p1, width=8, height=8)
  # END VIS CODE -------------------------------------------
} # end of levels
