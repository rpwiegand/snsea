library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(latex2exp)
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
  sigmaLabel  = paste("Gaussian Mutation $\\sigma = ", sigma, "$")
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


estimateConvergenceRate <- function(sigma=0.1, rhoMin=0.45, suffix="-nopop") {
  filename = paste("Results/unbounded-rv-s", sigma, "-r", rhoMin, suffix,".out", sep='') 
  
  print(paste("Reading file: ", filename), quote=F)
  results <- read.table(filename, header=T)
  aggregatedResults <- summarise(group_by(results, Generation), 
                                 AvgCoverPackingDiff=max(mean(CoverEpsilon) - mean(PackingEpsilon),0))
  aggregatedResults <- mutate(aggregatedResults,
                              ScaledCPDiff = AvgCoverPackingDiff/(max(AvgCoverPackingDiff)-min(AvgCoverPackingDiff)) - min(AvgCoverPackingDiff))
  return(sum(aggregatedResults$ScaledCPDiff))
}

estimateCoverPackingIntersection <- function(sigma=0.1, rhoMin=0.45, suffix="-nopop", maxHypotheticalGeneration=100000, minGen=10) {
  filename = paste("Results/unbounded-rv-s", sigma, "-r", rhoMin, suffix,".out", sep='') 
  print(paste("Reading file: ", filename), quote=F)
  results <- read.table(filename, header=T)
#  results <- filter(resultsRaw,
#                    CoverEpsilon < 3,
#                    PackingEpsilon > 0)
  
  
  coverModel <- loess(CoverEpsilon ~ Generation, results)
  hypotheticalCover <- predict(coverModel, newdata=data.frame(Generation=minGen:maxHypotheticalGeneration))

  packingModel <- loess(PackingEpsilon ~ Generation, results)
  hypotheticalPacking <- predict(packingModel, newdata=data.frame(Generation=minGen:maxHypotheticalGeneration))
  
  # LOESS is fine, if we do not have to extrapolate.  Otherwse we'll have to assume a specific
  # function.  If the curves have not intersected, then use log(x) for cover and sqrt(x) for packing
  trueMaxGen <- max(results$Generation)
  if (hypotheticalCover[trueMaxGen-minGen] > hypotheticalPacking[trueMaxGen-minGen]) {
    coverModel <- lm(CoverEpsilon ~ 1/Generation, results)
    hypotheticalCover <- predict(coverModel, newdata=data.frame(Generation=minGen:maxHypotheticalGeneration))
    
    packingModel <- lm(PackingEpsilon ~ sqrt(Generation), results)
    hypotheticalPacking <- predict(packingModel, newdata=data.frame(Generation=minGen:maxHypotheticalGeneration))
  }

  return(which(hypotheticalCover - hypotheticalPacking < 0)[[1]] + minGen)
}


plotUnboundedCoverPacking <- function(rhoMin=0.45, sigma=0.1, suffix="-mu2-lambda8", savePDF=False) {
  filename = paste("Results/unbounded-rv-s", sigma, "-r", rhoMin, suffix,".out", sep='') 
  #filename = paste("Results/unbounded-sig-", sigma, "-k-3-rho-", rhoMin, suffix,".out", sep='') 
  print(paste("Reading file: ", filename), quote=F)
  
  results <- read.table(filename, header=T)
  maxMinSpar <- max(results$MinArchiveSparseness)
  #minY <- min(results$MinArchiveSparseness)
  #maxY <- max(results$CoverEpsilon)
  
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
  minY = min(longAggregatedResults$value)
  maxY = max(longAggregatedResults$value)
  
  maxGen <- max(longAggregatedResults$Generation)
  rhoMinLabel = paste("$\\rho_{min} =", rhoMin, "$")
  sigmaLabel  = paste("Guassian Mutation $\\sigma = ", sigma, "$")
  titleLabel = paste("Unbounded 5D Euclidean Space")
  keyPoint <- findKeyPoint(results)
  
  p1 <- ggplot(longAggregatedResults, aes(x=Generation, y=value, group=variable, color=variable)) + 
    geom_hline(yintercept=rhoMin, size=0.5, color="black", linetype="longdash") +
    geom_line(size=1.25) +
    ylim(c(minY, maxY)) +             #+0.25
    annotate("text", x=-20, y=rhoMin+0.05, label=TeX(rhoMinLabel), hjust=0, size=4) +
    annotate("text", x=120, y=1.3, label="Packing continues to increase\nbut space is no longer\nbeing covered",
             vjust=1, hjust=0, size=3, color="darkblue") +
    #annotate("text", x=maxGen, y=1.6, label=TeX("Eventually converges to $\\epsilon$-Net"),
    #         vjust=1, hjust=1, size=4, color="firebrick") +
    scale_color_brewer(name="Archive\nMeasure",
                       palette="Set1",
                       labels=c("Cover\nEpsilon\n",
                                "Packing\nEpsilon\n", 
                                "Minimum\nSparseness")) +
    xlab("Generation") +
    ylab("Measure") +
    ggtitle(titleLabel, subtitle=TeX(sigmaLabel)) +
    theme(text=element_text(family="Times", size=16),
          legend.position="bottom")
  
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


barplotUnboundedConvergencePointByEval <- function(filename='Results/UnboundedConvergencePoint.csv', savePDF=F) {
  ucp <- read.csv(filename, header=T, sep='|')
  ucp <- mutate(ucp,
                EvaluationCount=GenerationConverged*PopulationSize,
                UsesPopulation=factor( c("mulambda","nopop")[ (ExperimentalGroup == 'NoPop')+1 ] ))
  
  p <- ggplot(ucp, aes(x=reorder(ExperimentalGroup, EvaluationCount), y=EvaluationCount, fill=UsesPopulation)) + 
        geom_bar(stat="identity") +
        scale_fill_manual(name="",
                          breaks=c("nopop","mulambda"),
                          values=c("steelblue","firebrick"),
                          labels=c("No Population","Population")) +
   # labels=c(TeX("$(\\mu,\\lambda)$ Population","No Population"))) +
  coord_flip() +
        xlab(TeX("Experimental Group, $(\\mu,\\lambda)$")) +
        ylab("Number of Evaluations before\nPacking & Cover Converged") +
        ggtitle("How Do Populations Affect SNS Converence?") + 
        theme(text=element_text(family="Times", size=16))
  
  print(p)
  
  if (savePDF) {
    pdfFilename = gsub("Results", "Graphs", gsub("csv", "pdf", filename))
    print(paste("Saving to file", pdfFilename))
    ggsave(pdfFilename, width=7, height=4, units="in")
  }
}


visualizePopulations <- function(gens, tr=0, filename='Results/viz-archive-and-pop.csv'){
  vizRaw <- read.csv(filename, header=T)
  vizFiltered <- mutate(filter(vizRaw, trial==tr, generation %in% gens),
                        isarchive=(whichPop!="archive"))
  
  p<- ggplot(vizFiltered, aes(x=x, y=y, color=whichPop, size=isarchive)) + 
    geom_point() +
    scale_color_manual(values=c("darkgray","red", "blue"),
                       name="") +
    facet_grid(. ~ generation, scales="fixed") +
    #facet_grid(generation ~ ., scales="fixed") +
    #facet_wrap(generation ~ .) + 
    scale_size_discrete(name="Population", range=c(1.25,2.5)) +
    ylab("") + xlab("") + # xlab("Generation") +
    ggtitle("Populations Do Not Cover Like Archives (by gen.)") +
    theme_bw() +
    theme(text=element_text(family="Times", size=18),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position="none") +
    guides(size=FALSE)
  
  print(p)
}



barplotUnboundedConvergencePoint <- function(filename='Results/UnboundedConvergencePoint.csv', savePDF=F) {
  ucp <- read.csv(filename, header=T, sep='|')
  ucp <- mutate(ucp,
                UsesPopulation=factor( c("mulambda","nopop")[ (ExperimentalGroup == 'NoPop')+1 ] ))
  
  p <- ggplot(ucp, aes(x=reorder(ExperimentalGroup, Order), y=GenerationConverged, fill=UsesPopulation)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(name="",
                      breaks=c("nopop","mulambda"),
                      values=c("firebrick","steelblue"),
                      labels=c("No Population","Population")) +
    # labels=c(TeX("$(\\mu,\\lambda)$ Population","No Population"))) +
    coord_flip() +
    xlab(TeX("Experimental Group, $(\\mu,\\lambda)$")) +
    ylab("Generation at which\nPacking & Cover Converged") +
    ggtitle("How Do Populations Affect SNS Converence?") + 
    theme(text=element_text(family="Times", size=16), legend.position="none")
  
  print(p)
  
  if (savePDF) {
    pdfFilename = gsub("Results", "Graphs", gsub("csv", "pdf", filename))
    print(paste("Saving to file", pdfFilename))
    ggsave(pdfFilename, width=7, height=4, units="in")
  }
}


plotChildCoverRatio <- function(filename='Results/unbounded-rv-s0.1-r0.45-pop2c16.out') {
  coverData <- read.table(filename, header=T)
  coverData <- mutate(coverData, ChildCoverRatio = ChildCover/CoverEpsilon )
  aggCoverData <- summarize(group_by(coverData, Generation), 
                            MeanChildCover=mean(ChildCover),
                            MeanChildCoverRatio=mean(ChildCoverRatio))
  #p <- ggplot(aggCoverData, aes(x=Generation, y=MeanChildCoverRatio)) + 
  #  geom_line(size=1.25, color="firebrick") +
  p <- ggplot(coverData, aes(x=Generation, y=ChildCover)) +
    #geom_line(size=0.5, color="darkgray") +
    geom_smooth(lwd=1.25, color="firebrick") +
    ylab("Cover Epsilon") +
    ggtitle("Children Population Cover",
            subtitle=TeX("(2,16)-SNSEA, 5D, $\\sigma = 0.45$, $\\rho_{min} = 0.1$")) +
    theme(text=element_text(family="Times", size=18))
  print(p)
}