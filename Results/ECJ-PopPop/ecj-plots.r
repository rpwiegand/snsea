library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(latex2exp)

darkCBPallet <- c(rgb(0.2873239, 0.7098592, 0.3971831),
                  rgb(0.5464789, 0.3971831, 0.4507042),
                  rgb(221/255, 170/255, 51/255))

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


findKeyPoint <- function(results) {
  idx <- which(results$CoverEpsilon < results$PackingEpsilon)[1]
  return( c(results$Generation[idx],results$CoverEpsilon[idx]) )
}

# =============================================================================
#  Data Access Functions
# =============================================================================


# -----------------------------------------------------------------------------
# Pull all the files for convergence results on the Euclidean space together
getRVResults <- function(sVals = c('01','02', '03'), 
                         rVals=c('02','04','06'), 
                         basename="boundedrv-500",
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



# -----------------------------------------------------------------------------
# Pull all the files for convergence results on the Euclidean space together
# for different sized populations.
getConvResults <- function(sVal = '01', 
                           rVal = '02', 
                           popSizes=c(2,4,8,16,32),
                           basename="boundedrv-conv",
                           suffix=".YY") {
  df <- NULL 
  for (m in popSizes) {
    for (l in popSizes) {
      if (l > m) {
        filename <- paste(basename, '-s', sVal, '-r', rVal, '-mu', m, 'lam', l, suffix, sep='')
        cat(paste("Reading file", filename, '\n'))
        tmpdf <- read.table(filename, 
                            col.names=c("Trial", "Cover", "GenConv", "MaxGens", "ConvergeInd"),
                            flush=T)
      
        numRows <- dim(tmpdf)[1]
        tmpdf$mu <- rep(as.numeric(m), numRows)
        tmpdf$llam <- rep(as.numeric(l), numRows)
      
        if (is.null(df)) df <- tmpdf
        else df <- rbind(df, tmpdf)
      }
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
          ggtitle("Bounded Euclidean (5D)") +
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
plotRVConvBarPlot <- function(df, titlePrefix="Bounded") {
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
    #scale_fill_brewer(palette="Set2", name=TeX("$\\sigma$")) +
    scale_fill_manual(values=darkCBPallet, name=TeX("$\\sigma$")) +
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
plotSmallMult500rv <- function(df, titlePrefix="Bounded", ylim=c(0.5,1.25)) {
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
          strip.text.y = element_text(size=7),          
          panel.background = element_blank()) 
  
  print(p)
  return(p)
}


# -----------------------------------------------------------------------------
#  Produce a plot of the Euclidean space convergence curves in terms of Cover,
#  Packing, and MinSparsness for a particular sigma and rhomin.
plotRV500CoverPacking <- function(df, titlePrefix="Bounded (Archive Selection)", inSigma=0.1, inRhoMin=0.2, 
                                  packingAnnotation=T, doFilter=T) {
  if (doFilter)
    results = filter(df, sigma==inSigma, rho==inRhoMin)
  else
    results = df
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
    #scale_color_brewer(name="Archive\nMeasure",
    #                   palette="Set1",
    #                   labels=c("Cover\nEpsilon\n",
    #                            "Packing\nEpsilon\n", 
    #                            "Minimum\nSparseness")) +
    scale_color_manual(name="Archive\nMeasure",
                       values=darkCBPallet,
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


plotDotWhiskerConvergenceRates <- function(df) {
  df2 <- summarise(group_by(df, mu, llam), 
                   AverageGensToConvergence=mean(GenConv),
                   MedianGensToConvergence=median(GenConv),
                   LowerCI=getLowerCI(GenConv),
                   UpperCI=getUpperCI(GenConv))
  #df2 <- mutate(df2,
  #              MuLab=paste("mu=", as.character(mu), sep=''))
  
  p <- ggplot(df2, aes(y=AverageGensToConvergence, x=factor(llam))) +
         geom_linerange(aes(ymin=LowerCI, ymax=UpperCI), size=2, color="firebrick") + 
         geom_point(aes(y=MedianGensToConvergence), size=4, shape=17) +
         geom_point(size=5, shape=21, color="firebrick", fill="firebrick") +
         #scale_x_discrete(breaks=c(2,4,8,16,32),
         #                  labels=c("2", "4", "8", "16", "32"),
        #                limits=c(0,34)) + 
         facet_grid(. ~ mu) +
         xlab(TeX("$\\lambda$")) +
         ylab("Generations to Convergences") +
         ggtitle("Convergence Results for Different Population Sizes",
                 subtitle=TeX("$\\mu$")) +
         theme_bw() +
         theme(text=element_text(family="Times", size=12),
               axis.line = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               #panel.border = element_blank(),
               panel.background = element_blank())   
         
 print(p)
 
 return(p)
}


visualizePopulations <- function(gens, tr=0, suffix='-popsel', ncol=6, titleStr="Selecting from Population (by Gen)") {
  filename=paste('viz-archive-and-pop-s01-r0001-mu16lam32',suffix,'.csv',sep='')
  vizRaw <- read.csv(filename, header=T)
  vizFiltered <- mutate(filter(vizRaw, trial==tr, generation %in% gens),
                        isarchive=(whichPop!="archive"))
  
  p<- ggplot(vizFiltered, aes(x=x, y=y, color=whichPop, size=isarchive)) + 
    geom_point(size=1.2, alpha=0.75) +
    scale_color_manual(values=c("darkgray","firebrick", "darkgray"),
                       name="") +
    facet_wrap(. ~ generation, scales="fixed", ncol=ncol) +
    #facet_grid(generation ~ ., scales="fixed") +
    #facet_wrap(generation ~ .) + 
    scale_size_discrete(name="Population", range=c(1.25,2.5)) +
    ylab("") + xlab("") + # xlab("Generation") +
    xlim(c(0,1)) + ylim(c(0,1)) +
    ggtitle(titleStr) +
    theme_bw() +
    theme(text=element_text(family="Times", size=12),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size=10),
          legend.position="none") +
    guides(size=FALSE)
  
  print(p)
}


# =============================================================================
#  Output / Saving Functions
# =============================================================================

produceAllPlots <- function(whichPlots=1:6) {
  df500 <- getRVResults()
  dfconv <- getConvResults()
  
  if (1 %in% whichPlots) {
    p1 <- plotSmallMult500rv(df500)
    ggsave("bounded-500-sm-mu2lam8.pdf", width=13, height=7, units="cm")
  }
  
  if (2 %in% whichPlots) {
    p2 <- plotRV500CoverPacking(df500)
    ggsave("bounded-500-s01r02-mu2lam8.pdf", width=13, height=7, units="cm")
  }
  
  if (3 %in% whichPlots) {
    p3 <- plotDotWhiskerConvergenceRates(dfconv)
    ggsave("bounded-conv-s02-r02.pdf", width=13, height=7, units="cm")
  }
  
  if (4 %in% whichPlots) {
    p4 <- visualizePopulations(seq(from=0, to=110, by=10), 0, '-popsel', 12)
    ggsave("viz-archive-and-pop-s01-r0001-mu16lam32-popsel.pdf", width=13, height=4, units="cm")
  }

  if (5 %in% whichPlots) {
    p5 <- visualizePopulations(seq(from=0, to=110, by=10), 0, '-archsel', 12, "Selecting from Archive (by Gen)")
    ggsave("viz-archive-and-pop-s01-r0001-mu16lam32-archsell.pdf", width=13, height=4, units="cm")
  }

  if (6 %in% whichPlots) {
    df <- read.table('boundedrv-500-vizexamp-mu2lam8-full.XX', header=T)
    p6 <- plotRV500CoverPacking(df, doFilter=F)
    ggsave("boundedrv-500-vizexamp-mu2lam8-full.pdf", width=13, height=7, units="cm")
  }
  
}