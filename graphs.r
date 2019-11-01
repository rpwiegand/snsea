library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(latex2exp)

results <- read.table('Results/results-sig-0.3-k-3-rho-0.6.out', header=T)
aggregatedResults <- summarise(group_by(results, Generation), 
                               AvgCover=mean(CoverEpsilon),
                               AvgPacking=mean(PackingEpsilon),
                               #AvgArchiveSize=mean(ArchiveSize),
                               AvgMinSparseness=mean(MinArchiveSparseness))
aggregatedArchiveSize <- summarise(group_by(results, Generation),
                               MinArchiveSize=min(ArchiveSize),
                               AvgArchiveSize=mean(ArchiveSize),
                               MaxArchiveSize=max(ArchiveSize))

longAggregatedResults <- filter(melt(aggregatedResults, id.vars=c("Generation")),
                                Generation < 500)

p1 <-ggplot(longAggregatedResults, aes(x=Generation, y=value, group=variable, color=variable)) + 
  geom_hline(yintercept=0.6, size=0.5, color="black", linetype="longdash") +
  annotate("text", x=210, y=0.65, label=TeX("$\\rho = 0.6$")) +
  geom_line(size=1.25) +
  scale_color_brewer(name="Archive Measure",
                     palette="Set1",
                     labels=c("Cover Epsilon","Packing Epsilon", "Min. Sparseness")) +
  xlab("Generation") +
  ylab("Measure") +
  theme(text=element_text(family="Times", size=14))

print(p1)

p2 <-ggplot(filter(results, Generation<500), aes(x=Generation, y=ArchiveSize)) + 
  geom_point(size=0.5, color="darkgray") + 
  geom_smooth(size=1.25, method="loess", color="firebrick", se=TRUE) +
  xlab("Generation") +
  ylab("Archive Size") +
  theme(text=element_text(family="Times", size=14))

print(p2)
