library(ggplot2)
library(dplyr)
library(viridis)
library(gridExtra)

setwd('~/sal_enhancers/coverage_impact')

#=======================================================================================================================
#     DFE
#=======================================================================================================================

# dfe <- read.csv('binned_dfe_allregions.csv')
dfe_31 <- subset(read.csv('~/sal_enhancers/dfe/binned_gammadfe_allregions_nes.csv'), region=='0fold' & rep==0)
dfe_31$n <- '31'

dfe_7 <- read.csv('binned_dfe_ss7_0fold.csv')
dfe_7$n <- '7'

dfe <- rbind(dfe_31, dfe_7)

levels(dfe$bin) <- c("0 - 1", "1 - 10", "10 - 100", "> 100")
str(dfe)

dfe_plot <- ggplot(dfe, aes(x=bin, y=proportion, fill=n)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = "bottom", legend.key.size = unit(.28, "cm"),
        legend.margin = margin(t = -1, r = 0, b = 2, l = -10, unit = "mm"),
        legend.text = element_text(colour = "black",size=6),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        plot.title = element_text(size=9, vjust=-6, hjust=0.007, face='bold'),
        strip.background = element_blank(),
        plot.margin = margin(t=0, r=5.5, b=5.5, l=5.5, unit = "pt"))+
  xlab(expression(N[e] * s)) + ylab('proportion of variants')


png('0fold_dfe_cov_comp.png', width=3, height=3, units='in', res=320)

dfe_plot

dev.off()


