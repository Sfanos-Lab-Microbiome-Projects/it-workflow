rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library("pheatmap")
library(viridis)
# ---------------------------------------------------
# configure per user ---
#setwd("it-workflow/code")
setwd(paste0(path.expand("~"),"/Desktop/it-workflow/code"))

analysisdir = "../analysis/P20-heatmap-species"
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
pdfdir = paste(analysisdir, "/pdfs/", sep="")
unlink(pdfdir, recursive=TRUE)
dir.create(pdfdir)
# ---------------------------------------------------
tablefile = "../analysis/P07-filtered-sum-taxa-db4.0/final/species.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)

meta           = A[,1:2]
meta           = data.frame(meta)
rownames(meta) = apply( meta[ , c("Region","SampleID") ], 1, paste, collapse = "." )
rownames(A)    = apply( A[ , c("Region","SampleID") ], 1, paste, collapse = "." )
A              = A[1:nrow(A),3:ncol(A)]
A              = t(A)

# have A exclude the rows that have metadata (1 to lastmdex)
A           = as.matrix(A)

# format row names in A to be shorter --
rownames(A) = gsub("k__Bacteria; ", "", rownames(A))
rownames(A) = gsub(";", ".", rownames(A))
rownames(A) = gsub(" ", "", rownames(A))
rownames(A) = gsub("__", "-", rownames(A))
rownames(A) = gsub("k_Bacteria.", "", rownames(A))

for (j in 1:nrow(A)){
  rownames(A)[j] = paste(unlist(strsplit(rownames(A)[j], ".g-", fixed = TRUE))[-1], collapse=".")
}
vals = array(0, dim(A))
for (i in 1:nrow(A)){
  for (j in 1:ncol(A)){
    vals[i,j] = as.numeric(as.character(A[i,j]))
  }
}
colnames(vals) = rownames(meta)
rownames(vals) = rownames(A)

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

# color scheme for % values in the heatmap --------
colscheme <- c("#266491", "#1b90e3", "#5eb4f0", "#b1dbf9", "#ffffff")
CnormedCOLORS  = colorRampPalette(colscheme)(20)
CnormedCOLORS2 = c(viridis(20), viridis(20)[20], viridis(20)[20],viridis(20)[20])
colscheme <- c("#f3d731", "#ebe798", "#44b2d4", "#036a8a")
CnormedCOLORS3  = rev(colorRampPalette(colscheme)(20))
colscheme <- c("#ffffff", "#98d7eb", "#44b2d4", "#036a8a")
CnormedCOLORS4  = colorRampPalette(colscheme)(20)
CnormedCOLORS5  = colorRampPalette(rev(brewer.pal(20,"YlGnBu")))(20)

# adding a pseudocount to improve log-scale visualization
vals1 = log10((vals/100)+0.001)

pdf(file=paste(pdfdir, "ordered.1.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "clustered.1.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "ordered.2.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS2, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "clustered.2.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS2, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()


pdf(file=paste(pdfdir, "ordered.3.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS3, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "clustered.3.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS3, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "ordered.4.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS4, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "clustered.4.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS4, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "ordered.5.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS5, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()

pdf(file=paste(pdfdir, "clustered.5.pdf", sep=""), width=7, height=5)
pheatmap(vals1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=meta, fontsize_row=4, fontsize = 6, annotation_colors=mycolors, scale="none", color=CnormedCOLORS5, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()
