library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(latex2exp)

darkCBPallet <- c(rgb(0.2873239, 0.7098592, 0.3971831),
                  rgb(0.5464789, 0.3971831, 0.4507042),
                  rgb(221/255, 170/255, 51/255))




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





compileArchive <- function(gen, trials) {
  runningFilename <- paste("vizarchive-", gen, "-s01-r02-mu1lam4.XX", sep='')
  rtDF <- read.table(runningFilename, header=T)
  
  #vizarchive-50-s01-r02-mu1lam4.XX
  baseFilename <- paste("archive", gen, "/archive",gen, sep='')
  
  archive <- NULL
  for (trial in trials) {
    filename <- paste(baseFilename, trial, ".csv", sep='')
    tmpArch <- read.csv(filename, header=F)
    colnames(tmpArch) <- c("x", "y", "Index")
    
    n <- dim(tmpArch)[1]

    # Grab the cover epislon estimate for the last generation and the half-way point    
    # Use this to estimate whether or not the archive likely converged, and that as a column
    lastCover <- filter(rtDF, Trial==trial, Generation==(gen-1))$CoverEpsilon
    halfCover <- filter(rtDF, Trial==trial, Generation==floor(gen/2))$CoverEpsilon
    if ( (lastCover - halfCover) < 0.00001) tmpArch$Converged <- factor(rep('Y', n), levels=c('Y','N'))
    else                                    tmpArch$Converged <- factor(rep('N', n), levels=c('Y','N'))
    
    # Add columns for the trial and the generation we ran to
    tmpArch$Trial <- rep(trial, n)
    tmpArch$Generation <- rep(gen, n)
    
    # Grow the dataset
    if (is.null(archive)) archive <- tmpArch
    else archive <- rbind(archive, tmpArch)
  }
  
  return(archive)
}


assembleArchives <- function(genCounts, trials) {
  #x <- NULL
  #y <- NULL
  #Generation <- NULL
  
  archive <- NULL
  for (gen in genCounts) {
    tmpArch <- compileArchive(gen, trials)
    if (is.null(archive)) archive <- tmpArch
    else archive <- rbind(archive, tmpArch)
    
    #filename <- paste("archive", gen, "/archive", gen, "-all.csv", sep='')
    #archive <- read.csv(filename, header=F)
    #x <- c(x, archive$V1)
    #y <- c(y, archive$V2)
    #Generation <- c(Generation, rep(gen, length(archive$V1)))
  }
  
  #return(data.frame(x,y,Generation))
  return(archive)
}


vizArchiveFill <- function(archive, maxLength=11, filterConverged=T) {
  r <- 0.6 #(maxLength-1) / maxLength
  theta <- seq(from=0, to=2*pi, length=100)
  circleDF <- data.frame(x=r*cos(theta),y=r*sin(theta))
  
  p <- NULL
  panelLabelObject <- element_text(size = 18)
  
  # Only use the archives that did *not* converge
  if (filterConverged) {
    archiveFiltered <- filter(archive, Converged == 'N')
    p <- ggplot(archiveFiltered, aes(x=x, y=y, color=factor(Trial))) +
           #geom_density_2d(aes(alpha = -(..level..)), bins=6, alpha=0.75) +
           geom_point(size=0.5, alpha=0.35) +
           #stat_density_2d(geom = "polygon", aes(alpha = (..level..), color=factor(Trial) ) ) + 
           ggtitle("Archive Behavior Points for Different Generation Converged Runs") 
    panelLabelObject <- element_blank()
  }
  
  # Use all the archives!
  else {
    archiveFiltered <- archive
    p <- ggplot(archiveFiltered, aes(x=x, y=y)) +
           geom_point(size=0.5, alpha=0.35, color="firebrick") +
           stat_density_2d(geom = "polygon", aes(alpha = (..level..) ), bins=4) + 
           ggtitle("Archive Behavior Points for Different Generation All Runs") 
  }
  
  p <- p + facet_grid(. ~ Generation) +
    #geom_density_2d(aes(alpha = -(..level..)), bins=5) + #alpha=0.75) +
    #geom_point(size=0.5, alpha=0.35) + #, color="firebrick") + 
    #stat_density_2d(geom = "polygon", aes(alpha = (..level..) ), bins=4) + 
    scale_alpha_continuous(range = c(0, 1)) +    
    geom_path(data=circleDF, aes(x=x,y=y), size=1.0, color="black") +
    xlab("") + ylab("") + 
    theme_bw() +
    theme(text=element_text(family="Times", size=14),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.text.x = panelLabelObject,
          legend.position="none") 
  
  print(p)
}




# =============================================================================
#  Output / Saving Functions
# =============================================================================

produceAllPlots <- function(whichPlots=1:2) {
  
  df <- assembleArchives(c(50,100,200), 0:69)
  
  if (1 %in% whichPlots) {
    p1 <- vizArchiveFill(df, filterConverged = T)
    ggsave("vizarchive-conv-s01-r02-mu1lam4.pdf", width=21, height=8, units="cm")
  }
  
  if (2 %in% whichPlots) {
    p2 <- vizArchiveFill(df, filterConverged = F)
    ggsave("vizarchive-all-0s01-r02-mu1lam4.pdf", width=21, height=9, units="cm")
  }
  
}