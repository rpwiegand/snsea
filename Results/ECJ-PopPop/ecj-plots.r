library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(latex2exp)

# =============================================================================
#  Convenience Utility Functions
# =============================================================================

# -----------------------------------------------------------------------------
# Convenience functions to (safely) get the lower/upper CI for drawing
# error bounds given a vector.  If all values are the same, these return
# the mean.
getLowerCI <- function(vct) {
  if (length(unique(vct)) > 1)
    return (t.test(vct)$conf[1])
  else
    return (mean(vct))
}

getUpperCI <- function(vct) {
  if (length(unique(vct)) > 1)
    return (t.test(vct)$conf[2])
  else
    return (mean(vct))
}


# =============================================================================
#  Data Access Functions
# =============================================================================


# -----------------------------------------------------------------------------
# Pull all the files for convergence results on the Euclidean space together
getRVResults <- function(sVals = c('01','02', '03'), 
                         rVals=c('02','04','06'), 
                         basename="unboundedrv-500",
                         popSuffix="-mu2lam8",
                         suffix=".XX") {
  df <- NULL 
  for (s in sVals) {
    for (r in rVals) {
      filename <- paste(basename, '-s', s, '-r', r, popSuffix, suffix, sep='')
      cat(paste("Reading file", filename, '\n'))
      tmpdf <- read.table(filename, header=T, flush=T)
      
      numRows <- dim(tmpdf)[1]
      tmpdf$sigma <- rep(as.numeric(s)/10, numRows)
      tmpdf$rho <- rep(as.numeric(r)/10, numRows)
      
      if (is.null(df)) df <- tmpdf
      else df <- rbind(df, tmpdf)
    }
  }
  
  return(df)
}



# =============================================================================
#  Graphing Functions
# =============================================================================


# -----------------------------------------------------------------------------
# Create a heatmap showing ho many trials converged for different 
# parameterization.
plotBoundedConvHeatMap <- function(df) {
  aggDF <- summarise(group_by(df, sigma, rho), ConvergeCount=table(ConvergeFlag)[1])
  
  p <- ggplot(aggDF, aes(x=sigma, y=rho, fill=ConvergeCount)) + #, color="white")) +
          geom_tile(color="white", size=0.75) +
          xlab(TeX("$\\sigma$")) +
          ylab(TeX("$\\rho_{min}$")) +
          scale_fill_continuous(name="Trials (out of 50)\nthat Converged") +
          ggtitle("Unbounded Euclidean (5D)") +
          theme_bw() +
          theme(text=element_text(family="Times", size=14),
                axis.line = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()) 

  print(p)
  return(p)
}


# -----------------------------------------------------------------------------
#  Produce a grouped bar plot of the real value parameterizations showing
#  the percentagle of trials converged.
plotRVConvBarPlot <- function(df, titlePrefix="Unbounded") {
  aggDF <- summarise(group_by(df, sigma, rho), ConvergeCount=table(ConvergeFlag)[1])
  
  p <- ggplot(aggDF, aes(x=rho, y=100*ConvergeCount/50, fill=factor(sigma))) + #, color="white")) +
    geom_bar(stat="identity", position="dodge", color="white") +
    geom_hline(yintercept=seq(from=0, to=100, by=25), size=0.75, color="white") +
    geom_hline(yintercept=seq(from=0, to=100, by=5), size=0.25, color="white") +
    #xlab(TeX("$\\sigma$")) +
    xlab(TeX("$\\rho_{min}$")) +
    ylab("Percente of Trials Converged") +
    #scale_fill_continuous(name="Trials (out of 50)\nthat Converged") +
    #scale_fill_brewer(palette="Set1", name=TeX("$\\rho_{min}$")) +
    scale_fill_brewer(palette="Set2", name=TeX("$\\sigma$")) +
    ggtitle(paste(titlePrefix,"Euclidean (5D) Convergence ")) +
    theme_bw() +
    theme(text=element_text(family="Times", size=12),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 

  print(p)
  return(p)
}



# -----------------------------------------------------------------------------
#  Produce a plot of the Euclidean space convergence curves in terms of Cover
#  as a small multiple of line plots over the different sigma and rho values.
plotSmallMult500rv <- function(df, titlePrefix="Unbounded", ylim=c(0.5,1.25)) {
  # New facet label names for sigma variable
  sigma.labs <- paste("sigma", c(0.1, 0.2, 0.3), sep='=')
  names(sigma.labs) <- c(0.1,0.2,0.3)
  
  # New facet label names for rhomin variable
  rho.labs <- paste("rho_min", c(0.2,0.4,0.6), sep="=")
  names(rho.labs) <- c(0.2,0.4,0.6)  
  
  # Aggregated dataset
  smDF <- summarise(group_by(df, Generation, sigma, rho), 
                    AverageCover=mean(CoverEpsilon),
                    LowerCICover=getLowerCI(CoverEpsilon),
                    UpperCICover=getUpperCI(CoverEpsilon))
  
  # Make the actual plot
  p <- ggplot(smDF, aes(x=Generation, y=AverageCover)) + 
    scale_x_continuous(breaks=c(0,250,500), labels=c("","250", "500")) +
    facet_grid( rho ~ sigma,
                labeller = labeller(sigma=sigma.labs, rho=rho.labs)) +
    #geom_segment(aes(x=Generation, xend=Generation, y=LowerCICover, yend=UpperCICover), 
    #             lineend="round",
    #             color="gray") + 
    geom_line(color="firebrick", size=0.75) + 
    ggtitle(paste(titlePrefix, "5D Euclidean Cover Epsilon")) +
    ylab("Average Cover Epsilon") +
    #ylim(ylim) +
    theme_bw() +
    theme(text=element_text(family="Times", size=12),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  print(p)
  return(p)
}


# -----------------------------------------------------------------------------
#  Produce a plot of the Euclidean space convergence curves in terms of Cover,
#  Packing, and MinSparsness for a particular sigma and rhomin.
plotRV500CoverPacking <- function(df, titlePrefix="Unbounded", inSigma=0.1, inRhoMin=0.2, packingAnnotation=T) {
  results = filter(df, sigma==inSigma, rho==inRhoMin)
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
  minY = min(longAggregatedResults$value)
  maxY = max(longAggregatedResults$value)
  
  maxGen <- max(longAggregatedResults$Generation)
  rhoMinLabel = paste("$\\rho_{min} =", inRhoMin, "$")
  sigmaLabel  = paste("Guassian Mutation $\\sigma = ", inSigma, "$")
  titleLabel = paste(titlePrefix, "5D Euclidean Space")
  keyPoint <- findKeyPoint(results)
  
  p1 <- ggplot(longAggregatedResults, aes(x=Generation, y=value, group=variable, color=variable)) + 
    geom_hline(yintercept=inRhoMin, size=0.5, color="black", linetype="longdash") +
    geom_line(size=1.25) +
    ylim(c(minY, maxY)) +
    annotate("text", label=TeX(rhoMinLabel), size=3,
             x=500, y=inRhoMin+0.05, hjust=1, vjust=0) +
    scale_color_brewer(name="Archive\nMeasure",
                       palette="Set2",
                       labels=c("Cover\nEpsilon\n",
                                "Packing\nEpsilon\n", 
                                "Minimum\nSparseness")) +
    xlab("Generation") +
    ylab("") +
    ggtitle(titleLabel, subtitle=TeX(sigmaLabel)) +
    theme_bw() + 
    theme(text=element_text(family="Times", size=12))

  if (packingAnnotation) {
    p1 <- p1 + annotate("text", x=100, y=3.4,
                        label="Packing continues to increase\nbut space is no longer\nbeing covered",
                        vjust=1, hjust=0, size=2, color="darkorange")
  }
  
  print(p1)
  return(p1)  
}


polynomialGenerator <- function(coefficients, tVals=seq(from=-pi,to=pi,length=20)) {
  yVals <- NULL
  n = length(coefficients) - 1
  
  for (theta in tVals) {
    y <- 0
    for (idx in 0:n) {
      if (idx %% 2) ct <- cos(theta)
      else          ct <- sin(theta)
      y <- y + coefficients[idx+1] * ct^idx 
    }
    yVals <- c(yVals, y)
  }
  
  return(yVals)
}


# =============================================================================
#  Output / Saving Functions
# =============================================================================

produceAllPlots <- function(whichPlots=1:9) {
  if (1 %in% whichPlots) {
    p1 <- plotRVConvBarPlotNOPOP(getRVResults(basename="boundedrv-conv"), "Bounded")
    ggsave("bounded-conv-NOPOP.pdf", width=13, height=10, units="cm")
  }
  
  if (2 %in% whichPlots) {
    p2 <- plotRVConvBarPlotNOPOP(getRVResults(basename="unboundedrv-conv"), "Unbounded")
    ggsave("unbounded-conv-NOPOP.pdf", width=13, height=10, units="cm")
  }
  
  if (3 %in% whichPlots) {
    p3 <- plotRV500CoverPacking(getRVResults(basename="boundedrv-500", suffix=".XX"), 
                                titlePrefix="Bounded", inSigma=0.1, inRhoMin=0.2, packingAnnotation=F)
    ggsave("bounded-s01-r02-NOPOP.pdf", width=13, height=7, units="cm")
  }
  
  if (4 %in% whichPlots) {
    p4 <- plotRV500CoverPacking(getRVResults(basename="unboundedrv-500", suffix=".XX"), 
                                titlePrefix="Unbounded", inSigma=0.1, inRhoMin=0.2, packingAnnotation=F)
    ggsave("unbounded-s01-r02-NOPOP.pdf", width=13, height=7, units="cm")
  }
  
  if (5 %in% whichPlots) {
    p5 <- plotRV500CoverPacking(getRVResults(basename="unboundedrv-500", suffix=".XX"), 
                                titlePrefix="Unbounded", inSigma=0.2, inRhoMin=0.6, packingAnnotation=T)
    ggsave("unbounded-s02-r06-NOPOP.pdf", width=13, height=7, units="cm")
  }
  
  if (6 %in% whichPlots) {
    p6 <- plotSmallMult500rvNOPOP(getRVResults(basename="boundedrv-500", suffix=".XX"), "Bounded")
    ggsave("bounded-500sm-NOPOP.pdf", width=13, height=13, units="cm")
  }
  
  if (7 %in% whichPlots) {
    p7 <- plotSmallMult500rvNOPOP(getRVResults(basename="unboundedrv-500", suffix=".XX"), "Unbounded")
    ggsave("unbounded-500sm-NOPOP.pdf", width=13, height=13, units="cm")
  }
  
  if (8 %in% whichPlots) {
    p8 <- plotSmallMult500HammingNOPOP(getHammingResults(basename="hamming-500", suffix=".XX"))
    ggsave("hamming-500sm-NOPOP.pdf", width=13, height=5, units="cm")
  }
  
  if (9 %in% whichPlots) {
    plotHamming500CoverPacking(getHammingResults(basename="hamming-500", suffix=".XX"))
    ggsave("hamming-500-n10-NOPOP.pdf", width=13, height=7, units="cm")
  }
  
}