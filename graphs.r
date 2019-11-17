library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(latex2exp)

findKeyPoint <- function(results) {
  idx <- which(results$CoverEpsilon < results$PackingEpsilon)[1]
  return( c(results$Generation[idx],results$CoverEpsilon[idx]) )
}

plotBoundedCoverPacking <- function(rhoMin=0.6, sigma=0.3, savePDF=False) {
  filename=paste("Results/bounded-sig-", sigma, "-k-3-rho-", rhoMin, ".out", sep='') 

  results <- read.table(filename, header=T)
  maxMinSpar <- max(results$MinArchiveSparseness)
  
  aggregatedResults <- summarise(group_by(results, Generation), 
                                 AvgCover=mean(CoverEpsilon),
                                 AvgPacking=mean(PackingEpsilon),
                                 AvgMinSparseness=mean(MinArchiveSparseness))

  longAggregatedResults <- filter(melt(aggregatedResults, id.vars=c("Generation")),
                                  Generation < 500)
  longAggregatedResults <- mutate(longAggregatedResults,
                                  variable = factor(variable,
                                                    levels=c("AvgCover",
                                                             "AvgPacking",
                                                             "AvgMinSparseness"),
                                                    ordered=T))

  rhoMinLabel = paste("$\\rho_{min} =", rhoMin, "$")
  sigmaLabel  = paste("Guassian Mutation $\\sigma = ", sigma, "$")
  titleLabel = paste("Bounded 3D Euclidean Space")
  #keyPoint <- findKeyPoint(results)
  keyPoint <- c(25, 0.75)
  
  p1 <- ggplot(longAggregatedResults, aes(x=Generation, y=value, group=variable, color=variable)) + 
    geom_hline(yintercept=0.6, size=0.5, color="black", linetype="longdash") +
    geom_line(size=1.25) +
    geom_curve(x=keyPoint[1]+75, y=keyPoint[2]+0.5, 
               xend=keyPoint[1], yend=keyPoint[2], 
               arrow=arrow(length=unit(0.015,"npc")), 
               angle=270, 
               curvature=0.25,
               color="darkgray") +
    #
    annotate("text", x=500, y=0.65, label=TeX(rhoMinLabel), hjust=1, size=5) +
    annotate("text", keyPoint[1]+75, keyPoint[2]+0.5, 
             label=TeX("Archive Converged to $\\epsilon$-Net"), 
             hjust=0, color="black", size=4) +
    annotate("text", x=500, y=0.18, size=4, hjust=1,
             label="Nothing substantively\nnew is being added",
             color="firebrick") +
    #
    scale_color_brewer(name="Archive\nMeasure",
                       palette="Set1",
                       labels=c("Cover\nEpsilon\n",
                                "Packing\nEpsilon\n", 
                                "Minimum\nSparseness")) +
    xlab("Generation") +
    ylab("Measure") +
    ggtitle(titleLabel, subtitle=TeX(sigmaLabel)) +
    theme(text=element_text(family="Times", size=16))

  print(p1)
  
  if (savePDF) {
    pdfFilename = gsub("Results", "Graphs", gsub("out", "pdf", filename))
    print(paste("Saving to file", pdfFilename))
    ggsave(pdfFilename, width=5, height=5, units="in")
  }
}


plotUnboundedCoverPacking <- function(rhoMin=0.6, sigma=0.2, savePDF=False) {
  filename=paste("Results/unbounded-sig-", sigma, "-k-3-rho-", rhoMin, ".out", sep='') 
  
  results <- read.table(filename, header=T)
  maxMinSpar <- max(results$MinArchiveSparseness)
  minY <- min(results$MinArchiveSparseness)
  maxY <- max(results$CoverEpsilon)
  
  aggregatedResults <- summarise(group_by(results, Generation), 
                                 AvgCover=mean(CoverEpsilon),
                                 AvgPacking=mean(PackingEpsilon),
                                 AvgMinSparseness=mean(MinArchiveSparseness))
  
  longAggregatedResults <- filter(melt(aggregatedResults, id.vars=c("Generation")),
                                  Generation < 500)
  longAggregatedResults <- mutate(longAggregatedResults,
                                  variable = factor(variable,
                                                    levels=c("AvgCover",
                                                             "AvgPacking",
                                                             "AvgMinSparseness"),
                                                    ordered=T))
  
  maxGen <- max(longAggregatedResults$Generation)
  rhoMinLabel = paste("$\\rho_{min} =", rhoMin, "$")
  sigmaLabel  = paste("Guassian Mutation $\\sigma = ", sigma, "$")
  titleLabel = paste("Unounded 8D Euclidean Space")
  keyPoint <- findKeyPoint(results)
  
  p1 <- ggplot(longAggregatedResults, aes(x=Generation, y=value, group=variable, color=variable)) + 
    geom_hline(yintercept=0.6, size=0.5, color="black", linetype="longdash") +
    geom_line(size=1.25) +
    ylim(c(minY, maxY)) +
    annotate("text", x=maxGen, y=0.8, label=TeX(rhoMinLabel), hjust=1, size=5) +
    #annotate("text", x=170, y=1.9, label="Packing continues to increase\nbut space is no longer\nbeing covered",
    #         vjust=1, hjust=0, size=4, color="darkblue") +
    annotate("text", x=maxGen, y=1.8, label=TeX("Eventually converges to $\\epsilon$-Net"),
             vjust=1, hjust=1, size=4, color="firebrick") +
    scale_color_brewer(name="Archive\nMeasure",
                       palette="Set1",
                       labels=c("Cover\nEpsilon\n",
                                "Packing\nEpsilon\n", 
                                "Minimum\nSparseness")) +
    xlab("Generation") +
    ylab("Measure") +
    ggtitle(titleLabel, subtitle=TeX(sigmaLabel)) +
    theme(text=element_text(family="Times", size=16))
  
  print(p1)
  
  if (savePDF) {
    pdfFilename = gsub("Results", "Graphs", gsub("out", "pdf", filename))
    print(paste("Saving to file", pdfFilename))
    ggsave(pdfFilename, width=5, height=5, units="in")
  }
}



#p2 <-ggplot(filter(results, Generation<500), aes(x=Generation, y=ArchiveSize)) + 
#  geom_point(size=0.5, color="darkgray") + 
#  geom_smooth(size=1.25, method="loess", color="firebrick", se=TRUE) +
#  xlab("Generation") +
#  ylab("Archive Size") +
#  theme(text=element_text(family="Times", size=14))
#
#print(p2)
