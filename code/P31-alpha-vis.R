rm(list=ls()) #reset of all commands in system
library(RColorBrewer) #install these packages first
library(ggplot2)
library(grid)
library(reshape2)
# ---------------------------------------------------
# configure per user ---
#setwd("it-workflow/code") # git repository already cloned to desktop
setwd(paste0(path.expand("~"),"/Desktop/it-workflow/code"))

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


# updates -- 1. Re-work P31-alpha-vis to split samples by fresh vs frozen,
# compare, and add stats if there is a significant diff bw fresh vs frozen
# (only 3 samples per group)
meltA$Type = "fresh"
meltA[grepl("frozen", meltA$SampleID),"Type"] = "frozen"
meltA$PID = gsub("-(fresh|frozen)", "", as.character(meltA$SampleID))

mycolors$Type["fresh"]  = "#74c7b8"
mycolors$Type["frozen"] = "#749fc7"

outfile1 = paste(pdfdir, "alpha-diversity-by-type.2.pdf", sep="")
p1 <- ggplot(meltA, aes(x=Region, y=Value, shape=Type)) +
  geom_boxplot(mapping=aes(color=Region), alpha=1, outlier.size = NA, coef=1000) +
  geom_point(data=meltA, aes(fill=Region), color="black", position=position_dodge(width=0.75), alpha=0.65, size=2) +
  theme_bw() +
  scale_color_manual(values=mycolors$Region) +
  scale_fill_manual(values=mycolors$Region) +
  scale_shape_manual(values=c(21,23)) +
  theme(axis.text.x  = element_text(size=11, colour="black"),
        axis.text.y  = element_text(size=11, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.title.y = element_text(size=12, colour="black"),
        plot.title   = element_text(size=10, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio=0.85) +
  xlab(NULL) +
  ylab("Measure Value") +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  facet_wrap(~Measure, ncol=2, scales="free_y")
ggsave(outfile1, plot=p1, width=6, height=6)



# GLM analysis -----
for (measure in unique(meltA$Measure)){
  capture.output(print(paste0("Fixed Effects GLM results for ", measure)),
                file = paste(analysisdir,"/glm-stats-results.txt",sep=""), append=TRUE)
  capture.output(summary(glm(Value ~ Type + Region + PID, data = meltA[meltA$Measure == measure,])),
                file = paste(analysisdir,"/glm-stats-results.txt",sep=""), append=TRUE)
}
