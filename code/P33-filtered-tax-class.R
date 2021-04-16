rm(list=ls()) #reset of all commands in system
library(RColorBrewer) #install these packages first
library(ggplot2)
library(grid)
library(reshape2)
# ---------------------------------------------------
# configure per user ---
#setwd("it-workflow/code") # git repository already cloned to desktop
setwd(paste0(path.expand("~"),"/Desktop/it-workflow/code"))

analysisdir = "../analysis/P33-filtered-tax-class" #create the analysis directory you want outputs to go into
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
pdfdir = paste(analysisdir, "/pdfs/", sep="")
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)
# ---------------------------------------------------
tablefile = "../analysis/P17-sum-taxa-db4.0.pl/final/species.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
meltA <- melt(A, id.vars=c("Region", "SampleID"))
colnames(meltA) = c("Region", "SampleID", "Taxa", "PercentAbundance")
meltA$Type = "fresh"
meltA[grepl("frozen", meltA$SampleID),"Type"] = "frozen"
meltA$PID = gsub("-(fresh|frozen)", "", as.character(meltA$SampleID))

# select the 32 most abundant species ---
avgPerTaxa = aggregate(PercentAbundance ~Taxa, data=meltA, FUN=mean)
avgPerTaxa = avgPerTaxa[order(avgPerTaxa[,2], decreasing=TRUE)[1:32],]

meltA = meltA[meltA$Taxa %in% c(as.character(avgPerTaxa[,1])),]

meltA$Taxa = gsub("k__Bacteria; ", "", meltA$Taxa)
meltA$Taxa = gsub(";", ".",  meltA$Taxa)
meltA$Taxa = gsub(" ", "",   meltA$Taxa)
meltA$Taxa = gsub("__", "-", meltA$Taxa)
meltA$Taxa = gsub("k_Bacteria.", "", meltA$Taxa)

for (j in 1:nrow(meltA)){
  meltA[j,"Taxa"] = paste(unlist(strsplit(meltA[j,"Taxa"], ".g-", fixed = TRUE))[-1], collapse=".")
}

# add in "Other assignments"
for (myregion in unique(meltA$Region)){
  for (mysamp in unique(meltA$SampleID)){
    if (nrow(meltA[meltA$SampleID == mysamp & meltA$Region == myregion,])==0){
      next
    }
    otherAb = 100-sum(as.numeric(meltA[meltA$SampleID == mysamp & meltA$Region == myregion, "PercentAbundance"]))
    if (otherAb < 0){
      otherAb = 0
    }
    meltA = rbind(meltA, c(myregion, mysamp, "Other species", otherAb, meltA[meltA$SampleID == mysamp, "Type"][1], meltA[meltA$SampleID == mysamp, "PID"][1]))
  }
}

meltA$PercentAbundance  = as.numeric(as.character(meltA$PercentAbundance))

meltA$Taxa             = as.character(meltA$Taxa)
nonUnassigned          = unique(meltA$Taxa)[!unique(meltA$Taxa) %in% c("unassigned.s-unassigned", "Other species")]

avgPerTaxa = aggregate(PercentAbundance ~ Taxa, data=meltA[meltA$Taxa %in% nonUnassigned,], FUN=mean)
avgPerTaxa = avgPerTaxa[order(avgPerTaxa[,2], decreasing=TRUE),]
orderednonUnassigned   = as.character(avgPerTaxa[,1])

meltA$Taxa             = factor(meltA$Taxa, levels=rev(c(orderednonUnassigned, "unassigned.s-unassigned", "Other species")))

# color scheme --
GroupCols = c("#cbcbcb", "#eb3eab", "#FF9AA2", "#FFDAC1", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#2C3E50", "#2980B9", "#8E44AD", "#C0392B", "#F39C12", "#F1C40F", "#27AE60", "#196f3e")
numTx     = 34
GroupCols = rev(c(rev(colorRampPalette(GroupCols)(numTx))))

p1 <- ggplot(meltA, aes(x=SampleID, y=PercentAbundance, fill=Taxa)) +
geom_bar(position="stack", stat="identity", colour="white", size=0.25) +
theme_bw() +
scale_fill_manual(values=GroupCols) +
theme(axis.text.x = element_text(size=7, colour="black", angle=45, hjust=1, vjust=1),
      axis.text.y = element_text(size=10, colour="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text      = element_text(size=5),
      legend.title     = element_text(size=5),
      aspect.ratio     = 28) +
xlab(NULL)  +
ylab("% Abundance") +
facet_grid(~Region, scale="free_x", space="free_x")
outfile1 = paste(pdfdir, "percent-taxa-per-region-v-expected-per-sample.1.pdf", sep="")
ggsave(file=outfile1, plot=p1, width=10, height=5)
