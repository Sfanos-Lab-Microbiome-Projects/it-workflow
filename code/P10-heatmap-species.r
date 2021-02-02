rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library("pheatmap")
# ---------------------------------------------------
# configure per user ---
# setwd("it-workflow/code")
analysisdir = "../analysis/P10-heatmap-species"
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
pdfdir = paste(analysisdir, "/pdfs/", sep="")
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)
# ---------------------------------------------------
tablefile = "../analysis/P07-sum-taxa-db4.0.pl/final/species.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
# ---------------------------------------------------

# BEGIN color scheme for major sample metadata features of interest ---------
mycolors = c()
mycolors$Region["V2"]  = "#FF9AA2"
mycolors$Region["V3"]  = "#FFDAC1"
mycolors$Region["V4"]  = "#FFF7C1"
mycolors$Region["V67"] = "#E2F0CB"
mycolors$Region["V8"]  = "#B5EAD7"
mycolors$Region["V9"]  = "#C7CEEA"

mycolors$SampleID["atcc1"] = "#F2978B"
mycolors$SampleID["atcc2"] = "#D954A0"
mycolors$SampleID["atcc3"] = "#6D5FA6"
mycolors$SampleID["atcc4"] = "#518FBF"
mycolors$SampleID["atcc5"] = "#63BBBF"
# END color scheme for major sample metadata features of interest ---------

# color scheme for % values in the heatmap --------
colscheme <- c("gray92", "#2C3E50", "#2980B9", "#8E44AD", "#C0392B", "#D35400", "#F39C12", "#F1C40F", "#27AE60", "#196f3e")
CnormedCOLORS = colorRampPalette(colscheme[1:10])(20)[c(1,3:20)]

# selecting samples with metadata

pdf(file=paste(pdfdir, "ordered.pdf", sep=""), width=10, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "clustered.pdf", sep=""), width=10, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()
