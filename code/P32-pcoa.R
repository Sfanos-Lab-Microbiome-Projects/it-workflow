rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(grid)
library(MASS)
library(vegan)
library(ape)

# set the working directory to the exact code base location ---
#setwd("it-workflow/code")
setwd(paste0(path.expand("~"),"/Desktop/it-workflow/code"))

# set the location place for the analysis outputs ---
analysisdir = "../analysis/P32-pcoa/"
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)

# load meta-data ---
tablefile = "../analysis/P18-sum-alpha/alpha-diversity-per-sample-region.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)
meta = A[,c(1,2)]
meta$Type = "fresh"
meta[grepl("frozen", meta$SampleID),"Type"] = "frozen"
meta$PID = gsub("-(fresh|frozen)", "", as.character(meta$SampleID))

imeta = unique(meta[,c(2,3,4)])
rownames(imeta) = imeta[,1]

regions    = c("v2", "v3", "v4", "v67", "v8", "v9")
dist.types = c("bray-curtis", "jaccard", "unweighted-unifrac", "weighted-unifrac")
# loop through regions and distance types for a comprehensive analysis ---
for (region in regions){
  for (dist in dist.types){
    # specify distance file and output prefix ---
    distfile   = paste("../analysis/P15-clust-tree-cm-goods/", region,"/", "core-metrics-results/", dist, "-dm/distance-matrix.tsv", sep="")
    prefix     = paste(dist, "-", region, sep="")
    print(paste("Analyzing ", prefix, "...", sep=""))

    # load and refine distance matrix ---
    myDist           <- read.table(distfile, header=TRUE, sep="\t", check.names=FALSE)
    rownames(myDist) = myDist[,1]
    # limit distance matrix to IDs matching the sample.ids in meta ---
    myDist           = myDist[rownames(myDist) %in% imeta[,1],  colnames(myDist) %in% imeta[,1] ]
    # limit meta data frame to those IDs matching the remaining IDs in myDist
    thismeta         = imeta[imeta$SampleID %in% rownames(myDist),]

    # execute PCoA ---
    bc.pcoa   = pcoa(myDist)
    # the output in bc.pcoa$vectors has the coordinates in pcoa space ---
    # extract the main two principal coordinates from PCoA results ---
    myDistf     = data.frame(PC1 = bc.pcoa$vectors[,1], PC2 = bc.pcoa$vectors[,2])

    # extract variance components from PCoA results (for axes in the ggplot)---
    bc.Ax1PrctVar = sprintf("%3.2f",100*bc.pcoa$values$Relative_eig[1])
    bc.Ax2PrctVar = sprintf("%3.2f",100*bc.pcoa$values$Relative_eig[2])

    # match metadata to myDistf
    myDistf             = cbind(thismeta[match(rownames(myDistf), thismeta[,1]),], myDistf)
    rownames(myDistf)   = myDistf$SampleID

    # --------------------------------------------------------
    # BEGIN color scheme ---
    mycolors = c()

    mycolors$Type["fresh"]  = "#74c7b8"
    mycolors$Type["frozen"] = "#749fc7"
    # END color scheme ---
    # --------------------------------------------------------

    # PERMANOVA analysis -----
    capture.output(adonis(myDist ~ Type + PID, data = myDistf, permutations = 1000),
                         file = paste(analysisdir, prefix,".PERMANOVA-results.txt",sep=""))

    # color by sample name
    p1 <- ggplot(myDistf, aes(x=PC1,y=PC2)) +
    geom_path(aes(group=PID), color="#c8c8c8", alpha=0.25) +
    geom_point(aes(fill=Type), pch=21, color="black", size=4, alpha=0.7) +
    scale_fill_manual(values=mycolors$Type)  +
    xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
    ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    # save p1 as a pdf
    ggsave(paste(analysisdir, prefix,".pcoa.01.pdf",sep=""), p1, width=4, height=4)
  } # end of dist types loop
} # end of region loop

# ---------------------------------------------------
# begin multi-region PCoA analysis here
# ---------------------------------------------------
tablefile = "../analysis/P17-sum-taxa-db4.0.pl/final/species.txt"
A <- read.table(tablefile, sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)

meta           = A[,1:2]
meta           = data.frame(meta)
rownames(meta) = apply( meta[ , c("Region","SampleID") ], 1, paste, collapse = "." )

meta$Type = "fresh"
meta[grepl("frozen", meta$SampleID),"Type"] = "frozen"
meta$PID = gsub("-(fresh|frozen)", "", as.character(meta$SampleID))

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

myvegdists = c("bray", "jaccard", "gower", "canberra", "euclidean", "kulczynski")
for (mydist in myvegdists){
  # create distance matrix
  # between all samples across regions
  fullDist = vegdist(t(vals), method=,mydist, upper=TRUE)
  # execute PCoA ---
  bc.pcoa   = pcoa(fullDist)
  # the output in bc.pcoa$vectors has the coordinates in pcoa space ---
  # extract the main two principal coordinates from PCoA results ---
  myDistf     = data.frame(PC1 = bc.pcoa$vectors[,1], PC2 = bc.pcoa$vectors[,2])
  rownames(myDistf) = rownames(bc.pcoa$vectors)
  # extract variance components from PCoA results (for axes in the ggplot)---
  bc.Ax1PrctVar = sprintf("%3.2f",100*bc.pcoa$values$Relative_eig[1])
  bc.Ax2PrctVar = sprintf("%3.2f",100*bc.pcoa$values$Relative_eig[2])

  # match metadata to myDistf
  myDistf             = cbind(meta[match(rownames(myDistf), rownames(meta)),], myDistf)
  # --------------------------------------------------------
  # BEGIN color scheme ---
  mycolors = c()
  mycolors$Region["v2"]  = "#FF9AA2"
  mycolors$Region["v3"]  = "#FFDAC1"
  mycolors$Region["v4"]  = "#FFF7C1"
  mycolors$Region["v67"] = "#E2F0CB"
  mycolors$Region["v8"]  = "#B5EAD7"
  mycolors$Region["v9"]  = "#C7CEEA"

  mycolors$PID["Patient1"] = "#98ddca"
  mycolors$PID["Patient2"] = "#ffd3b4"
  mycolors$PID["Patient3"] = "#ffaaa7"


  # END color scheme ---
  # --------------------------------------------------------
  # color by sample name
  p1 <- ggplot(myDistf, aes(x=PC1,y=PC2)) +
  geom_point(aes(fill=Region, shape=Type), color="black", size=4.5, alpha=0.7) +
  scale_fill_manual(values=mycolors$Region)  +
  scale_color_manual(values=mycolors$Region) +
  scale_shape_manual(values=c(21,23)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, "/allregions.", mydist, ".pcoa.01.pdf",sep=""), p1, width=4, height=4)

  p1 <- ggplot(myDistf, aes(x=PC1,y=PC2)) +
  geom_point(aes(fill=Region, shape=Type), color="black", size=4.5, alpha=0.7) +
  scale_fill_manual(values=mycolors$Region)  +
  scale_color_manual(values=mycolors$Region) +
  scale_shape_manual(values=c(21,23)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(aspect.ratio=1) +
  facet_wrap(~PID, ncol=2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, "/allregions.", mydist, ".pcoa.02.pdf",sep=""), p1, width=5, height=5)


  p1 <- ggplot(myDistf, aes(x=PC1,y=PC2)) +
  geom_path(aes(group=PID), color="#c8c8c8", alpha=0.5) +
  geom_point(aes(fill=PID, shape=Type), color="black", size=4.5, alpha=0.7) +
  scale_fill_manual(values=mycolors$PID)  +
  scale_color_manual(values=mycolors$PID) +
  scale_shape_manual(values=c(21,23)) +
  xlab(paste("PCoA Axis 1 (", bc.Ax1PrctVar, "% of Variation)", sep="")) +
  ylab(paste("PCoA Axis 2 (", bc.Ax2PrctVar, "% of Variation)", sep="")) +
  theme_bw() +
  theme(aspect.ratio=1) +
  facet_wrap(~Region, ncol=3) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes=list(shape=21)))
  ggsave(paste(analysisdir, "/allregions.", mydist, ".pcoa.03.pdf",sep=""), p1, width=6, height=6)

} # end of vegdists
