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
    annotate("text", x=500, y=0.7, label=TeX(rhoMinLabel), hjust=1, size=5) +
    annotate("text", keyPoint[1]+75, keyPoint[2]+0.5, 
             label=TeX("Archive Converged to $\\epsilon$-Net"), 
             hjust=0, color="black", size=4) +
    annotate("text", x=500, y=0.16, size=4, hjust=1,
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
    ggsave(pdfFilename, width=5, height=4, units="in")
  }
}


plotUnboundedCoverPacking <- function(rhoMin=0.6, sigma=0.2, savePDF=False) {
  #filename="Results/unbounded-sig-0.1-k-3-rho-0.4-population.out"
  filename = paste("Results/unbounded-sig-", sigma, "-k-3-rho-", rhoMin, "-pop2c8.out", sep='') 
  
  results <- read.table(filename, header=T)
  maxMinSpar <- max(results$MinArchiveSparseness)
  #minY <- min(results$MinArchiveSparseness)
  #maxY <- max(results$CoverEpsilon)
  
  aggregatedResults <- summarise(group_by(results, Generation), 
                                 AvgCover=mean(CoverEpsilon),
                                 AvgPacking=mean(PackingEpsilon),
                                 AvgMinSparseness=mean(MinArchiveSparseness))
  
  longAggregatedResults <- filter(melt(aggregatedResults, id.vars=c("Generation")),
                                  Generation < 100)
  longAggregatedResults <- mutate(longAggregatedResults,
                                  variable = factor(variable,
                                                    levels=c("AvgCover",
                                                             "AvgPacking",
                                                             "AvgMinSparseness"),
                                                    ordered=T))
  minY = min(longAggregatedResults$value)
  maxY = max(longAggregatedResults$value)
  
  maxGen <- max(longAggregatedResults$Generation)
  rhoMinLabel = paste("$\\rho_{min} =", rhoMin, "$")
  sigmaLabel  = paste("Guassian Mutation $\\sigma = ", sigma, "$")
  titleLabel = paste("Unounded 3D Euclidean Space")
  keyPoint <- findKeyPoint(results)
  
  p1 <- ggplot(longAggregatedResults, aes(x=Generation, y=value, group=variable, color=variable)) + 
    geom_hline(yintercept=rhoMin, size=0.5, color="black", linetype="longdash") +
    geom_line(size=1.25) +
    ylim(c(minY, maxY)) +
    annotate("text", x=maxGen, y=rhoMin-0.01, label=TeX(rhoMinLabel), hjust=1, size=5) +
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
    ggsave(pdfFilename, width=5, height=4, units="in")
  }
}



plotHammingCoverPacking <- function(rhoMin=3, n=16, savePDF=False) {
  filename=paste("Results/hamming-r", rhoMin, "-n", n, ".out", sep='') 

  results <- read.table(filename, header=T)
  
  maxMinSpar <- max(results$MinArchiveSparseness)
  minY <- min(results$MinArchiveSparseness)
  maxY <- max(results$CoverEpsilon)

  aggregatedResults <- summarise(group_by(results, Generation), 
                                 AvgCover=mean(CoverEpsilon),
                                 AvgPacking=mean(PackingEpsilon),
                                 AvgMinSparseness=mean(MinArchiveSparseness))

  longAggregatedResults <- filter(melt(aggregatedResults, id.vars=c("Generation")),
                                  Generation < 200)
  longAggregatedResults <- mutate(longAggregatedResults,
                                  variable = factor(variable,
                                                    levels=c("AvgCover",
                                                             "AvgPacking",
                                                             "AvgMinSparseness"),
                                                    ordered=T))

  maxGen <- max(longAggregatedResults$Generation)
  rhoMinLabel = paste("$\\rho_{min} =", rhoMin, "$")
  nLabel  = paste("Bit-flip mutation using $1/n$")
  titleLabel = paste("Hamming Space, n=", n)
  keyPoint <- findKeyPoint(results)

  p1 <- ggplot(longAggregatedResults, aes(x=Generation, y=value, group=variable, color=variable)) + 
    geom_hline(yintercept=0.6, size=0.5, color="black", linetype="longdash") +
    geom_line(size=1.25) + 
    ylim(c(minY, maxY)) +
    annotate("text", x=maxGen, y=1, label=TeX(rhoMinLabel), hjust=1, size=5) +
    annotate("text", x=maxGen, y=9.75, label=TeX("Converges to $\\epsilon$-Net"),
             vjust=1, hjust=1, size=4, color="firebrick") +
    scale_color_brewer(name="Archive\nMeasure",
                       palette="Set1",
                       labels=c("Cover\nEpsilon\n",
                                "Packing\nEpsilon\n", 
                                "Minimum\nSparseness")) +
    xlab("Generation") +
    ylab("Measure") +
    ggtitle(titleLabel, subtitle=TeX(nLabel)) +
    theme(text=element_text(family="Times", size=16))

  print(p1)

  if (savePDF) {
    pdfFilename = gsub("Results", "Graphs", gsub("out", "pdf", filename))
    print(paste("Saving to file", pdfFilename))
    ggsave(pdfFilename, width=5, height=4, units="in")
  }
}
