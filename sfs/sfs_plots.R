library(ggplot2)
library(viridis)

setwd('~/sal_enhancers/sfs')

sfs_dat <- read.csv('plottable_sfs.csv')

raw <- subset(sfs_dat, bs_rep==0)

sfs <- ggplot(raw, aes(x=freq, y=proportion, fill=region)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_viridis(discrete = T) +
  xlab('Derived allele frequency')

png('sfs_plot.png', width=8, height=3, res=320, units='in')

sfs

dev.off()