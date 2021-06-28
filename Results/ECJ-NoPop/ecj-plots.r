library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(latex2exp)

getHammingConvResults <- function(nVals = c(6,8,10,12,14)) {
  basename = "hamming-conv-n"
  suffix = ".YY"
  
  df <- NULL 
  for (n in nVals) {
    filename <- paste(basename, n, suffix, sep='')
    cat(paste("Reading file", filename, '\n'))
    tmpdf <- read.table(filename, header=T)
    
    numRows <- dim(tmpdf)[1]
    tmpdf$n <- rep(n, numRows)
    
    if (is.null(df)) df <- tmpdf
    else df <- rbind(df, tmpdf)
  }
  
  return(df)
}


getRVConvResults <- function(sVals = c('01','02', '03'), rVals=c('02','04','06'), basename="boundedrvi-conv") {
  suffix = ".YY"
  
  df <- NULL 
  for (s in sVals) {
    for (r in rVals) {
      filename <- paste(basename, '-s', s, '-r', r, suffix, sep='')
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


getBoundedConvHeatMapNOPOP <- function() {
  df     <- getRVConvResults()
  aggDF <- summarise(group_by(foo, sigma, rho), ConvergeCount=table(ConvergeFlag)[1])
  
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
  #+
          #theme(text=element_text(family="times", size=14))
  print(p)
}
