library(ggplot2)
library(viridis)

setwd('~/sal_enhancers/sfs')

sfs_dat <- read.csv('plottable_sfs.csv')

raw <- subset(sfs_dat, bs_rep==0)

sfs <- ggplot(raw, aes(x=freq, y=proportion, fill=region)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_viridis(discrete = T, guide = guide_legend(nrow=2)) +
  xlab('derived allele count') +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = "bottom",
        legend.text = element_text(colour = "black",size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        strip.background = element_blank())

png('sfs_plot.png', width=8, height=4, res=320, units='in')

sfs

dev.off()