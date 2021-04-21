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

analysisdir = "../analysis/P34-diff-ab-fresh-frozen" # create the analysis directory you want outputs to go into
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ---------------------------------------------------
taxalevels = c("alpha-diversity",
               "phylum",
               "class",
               "order",
               "family",
               "genus",
               "species")

for (txlevel in taxalevels){
  # ---------------------------------------------------
  # create output dirs ---
  txdir = paste(analysisdir, "/", txlevel, sep="") #subdirectory in analysisdir
  unlink(txdir, recursive=TRUE)
  dir.create(txdir)
  # ---------------------------------------------------
  pdfdir = paste(txdir, "/pdfs/", sep="") #subdirectory in analysisdir
  unlink(pdfdir, recursive=TRUE)
  dir.create(pdfdir)
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

  # LIST CRITICAL COMPARISONS ---
  my_comparisons = list( c("fresh", "frozen") )
  # END CRITICAL COMPARISONS ---

  # create array of features to assess ---
  # colnames of Awmeta are too complicated at this point
  # e.g. "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales;
  # f__Pasteurellaceae; g__Avibacterium; s__Avibacterium_paragallinarum"
  # we will remove semicolons, spaces and simplify the names before visualization
  colnames(Awmeta)[3:ncol(Awmeta)] = gsub("k_Bacteria.", "", colnames(Awmeta)[3:ncol(Awmeta)])
  colnames(Awmeta)[3:ncol(Awmeta)] = gsub(";", ".", colnames(Awmeta)[3:ncol(Awmeta)])
  colnames(Awmeta)[3:ncol(Awmeta)] = gsub(" ", "", colnames(Awmeta)[3:ncol(Awmeta)])
  colnames(Awmeta)[3:ncol(Awmeta)] = gsub("__", ".", colnames(Awmeta)[3:ncol(Awmeta)])
  colnames(Awmeta)[3:ncol(Awmeta)] = gsub("exported-", "", colnames(Awmeta)[3:ncol(Awmeta)])
  colnames(Awmeta)[3:ncol(Awmeta)] = gsub("-vector", "", colnames(Awmeta)[3:ncol(Awmeta)])
  colnames(Awmeta)[3:ncol(Awmeta)] = gsub("-", ".", colnames(Awmeta)[3:ncol(Awmeta)])

  # completed modifying colnames of taxa ---
  taxaFeatures = colnames(Awmeta)[3:ncol(Awmeta)] # this is the list we will visualize

  # BEGIN STATS CODE -------------------------------------------
  # this will hold all stats results for Q1 / Q2
  txStatsResults = c()

  for (myFeature in taxaFeatures){
    if (grepl("(PID|Type)", myFeature)){
      next
    }
    print(paste("Evaluating feature: ", txlevel, " ", myFeature, sep=""))
    if (sum(Awmeta[,myFeature] > 0) < 10){
      next
    }

    # add pseudo count for visualization / stats analysis
    Awmeta[,myFeature] = Awmeta[,myFeature] + 0.01

    resultsRow = c(myFeature)

    for (region1 in levels(Awmeta$Region)){
      for (c1 in levels(Awmeta$Type)){
        nc1        = length(Awmeta[Awmeta$Region==region1 &
                                   Awmeta$Type==c1, myFeature])
        resultsRow = c(resultsRow, nc1)
      }
      for (c1 in levels(Awmeta$Type)){
        mnc1       = mean(Awmeta[Awmeta$Region==region1 &
                                 Awmeta$Type==c1, myFeature])
        resultsRow = c(resultsRow, mnc1)
      }
      for (compi in 1:length(my_comparisons)){
        mwP = t.test(log10(Awmeta[Awmeta$Region==region1 & Awmeta$Type==my_comparisons[[compi]][1], myFeature]),
                     log10(Awmeta[Awmeta$Region==region1 & Awmeta$Type==my_comparisons[[compi]][2], myFeature]), paired=TRUE)$p.value

        resultsRow = c(resultsRow, mwP)
      }
    } # end of regions

    # add glm results  --
    glm1 <- glm(as.formula(paste(myFeature, " ~ Type + Region + PID", sep="")), data=Awmeta,  family="gaussian")
    typeCoef   = tryCatch({summary(glm1)$coefficients[2,1]}, error = function(e) {"n/a"})
    typePval   = tryCatch({summary(glm1)$coefficients[2,4]}, error = function(e) {"n/a"})
    resultsRow = c(resultsRow, typeCoef, typePval)

    glm1 <- glm(as.formula(paste("log10(", myFeature, ") ~ Type + Region + PID", sep="")), data=Awmeta,  family="gaussian")
    typeCoef   = tryCatch({summary(glm1)$coefficients[2,1]}, error = function(e) {"n/a"})
    typePval   = tryCatch({summary(glm1)$coefficients[2,4]}, error = function(e) {"n/a"})
    resultsRow = c(resultsRow, typeCoef, typePval)

    txStatsResults = rbind(txStatsResults,resultsRow)

    # END STATS CODE -------------------------------------------

    # BEGIN VIS CODE -------------------------------------------
    # ggplot2 ---
    myTitle = paste(myFeature, "\n", sep="")

    p1<-c()
    if (grepl("alpha", txlevel)){
      give.n <- function(x,min){
        return(c(y=min*0.9, label = length(x)))
      }

      p1 <- ggplot(Awmeta, aes_string(x="Type", y=myFeature)) +
      geom_path(aes(group=PID), alpha=0.25, color="#c8c8c8") +
      geom_boxplot(aes(color=Type),fill=NA, outlier.size=0, coef=1e100, alpha=0.9) +
      geom_point(aes(color=Type), alpha=0.95, size=3, position=position_dodge(width=0.75)) +
      theme_bw() +
      scale_color_manual(values=mycolors$Type) +
      theme(axis.text.x  = element_text(size=8, colour="black", angle=45, hjust=1, vjust=1),
            axis.text.y  = element_text(size=8, colour="black"),
            axis.title.x = element_text(size=12, colour="black"),
            axis.title.y = element_text(size=12, colour="black"),
            plot.title   = element_text(size=9, colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1.7) +
      xlab(NULL) + # no label for x axis
      ylab(paste(txlevel, " measure", sep="")) +
      ggtitle(myTitle) +
      facet_wrap(~Region, ncol=5) +
      stat_summary(aes(color=Type), size=3, fun.data = give.n, geom = "text", fun.args = list(min = min(Awmeta[,myFeature])))
    }else{
      give.n <- function(x,min){
        return(c(y=log10(min*0.5), label = length(x)))
      }

      p1 <- ggplot(Awmeta, aes_string(x="Type", y=myFeature)) +
      geom_path(aes(group=PID), alpha=0.25, color="#c8c8c8") +
      geom_boxplot(aes(color=Type),fill=NA, outlier.size=0, coef=1e100, alpha=0.9) +
      geom_point(aes(color=Type), alpha=0.95, size=3, position=position_dodge(width=0.75)) +
      theme_bw() +
      scale_color_manual(values=mycolors$Type) +
      theme(axis.text.x  = element_text(size=8, colour="black", angle=45, hjust=1, vjust=1),
            axis.text.y  = element_text(size=8, colour="black"),
            axis.title.x = element_text(size=12, colour="black"),
            axis.title.y = element_text(size=12, colour="black"),
            plot.title   = element_text(size=9, colour="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1.7) +
      xlab(NULL) + # no label for x axis
      ylab("% Ab.") +
      ggtitle(myTitle) +
      facet_wrap(~Region, ncol=5) +
      scale_y_log10(breaks=c(0.01, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100)) +
      stat_summary(aes(color=Type), size=3, fun.data = give.n, geom = "text", fun.args = list(min = min(Awmeta[,myFeature])))
    }
    myFeatureName = substr(myFeature,0,200)
    ggsave(paste(pdfdir, myFeatureName, ".01.pdf", sep=""), plot=p1, width=7, height=4)
    # END VIS CODE -------------------------------------------
  } # end of features loop

  # finalize colnames for txStatsResults ---
  txStatsResultsColNames = c("Feature")
  for (r1 in levels(Awmeta$Region)){
    for (c1 in unique(Awmeta$Type)){
      txStatsResultsColNames = c(txStatsResultsColNames,
                                 paste(r1, ".", c1,".N",sep=""))
    }
    for (c1 in unique(Awmeta$Type)){
      txStatsResultsColNames = c(txStatsResultsColNames,
                                 paste(r1, ".", c1,".mean",sep=""))
    }
    for (compi in 1:length(my_comparisons)){
      txStatsResultsColNames = c(txStatsResultsColNames,
        paste(r1, ".tt.log.Pval.",my_comparisons[[compi]][1],".vs.",my_comparisons[[compi]][2],sep=""))
    }
  }
  txStatsResultsColNames   = c(txStatsResultsColNames, "glm.std.Type.Coef", "glm.std.Type.Pval", "glm.log.Type.Coef", "glm.log.Type.Pval")
  colnames(txStatsResults) = txStatsResultsColNames
  txStatsResults = data.frame(txStatsResults)
  # get FDR-adjusted Pvalues
  txStatsResults$glm.std.Type.adjPval = p.adjust(as.numeric(as.character(txStatsResults$glm.std.Type.Pval)), method="fdr")
  txStatsResults$glm.log.Type.adjPval = p.adjust(as.numeric(as.character(txStatsResults$glm.log.Type.Pval)), method="fdr")

  txStatsResults = txStatsResults[,c(1:21,22:23,26,24:25,27)]

  write.table(x=txStatsResults, file=paste(txdir,"/", txlevel, "-stats-results.csv", sep=""), col.names=TRUE, row.names=FALSE, sep=",")
} # end of files loop
