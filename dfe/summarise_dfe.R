library(ggplot2)
library(dplyr)
library(viridis)

setwd('/Users/henryjuho/sal_enhancers/dfe')

# dfe <- read.csv('binned_dfe_allregions.csv')
dfe <- read.csv('binned_gammadfe_allregions_nes.csv')

str(dfe)
dfe$region <- factor(dfe$region, levels=c('intergenic', 'intergenic-enhancers', 'intron', 'intron-enhancers', 'utr', 'utr-enhancers', 'cds', 'cds-enhancers', 'all-enhancers'))
levels(dfe$region) <- c('intergenic', 'H3K27ac peaks (intergenic)', 'intron', 'H3K27ac peaks (intron)', 'UTR', 'H3K27ac peaks (UTR)', 'CDS', 'H3K27ac peaks (CDS)', 'H3K27ac peaks (all)')

levels(dfe$bin) <- c("0 - 1", "1 - 10", "10 - 100", "> 100")
str(dfe)

raw <- subset(dfe, rep==0)
bs <- subset(dfe, rep!=0)

cis <- summarise(group_by(bs, region, bin),
                 lwr=quantile(proportion, 0.025), upr=quantile(proportion, 0.975))

plot_data <- full_join(raw, cis)

dfe_plot <- ggplot(plot_data, aes(x=bin, y=proportion, fill=region)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_manual(values=viridis(10)[2:10]) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.5) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = "bottom", legend.key.size = unit(.28, "cm"),
        legend.margin = margin(t = -1, r = 0, b = 2, l = -10, unit = "mm"),
        legend.text = element_text(colour = "black",size=6),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        strip.background = element_blank())+
  xlab(expression(N[e] * s)) +
  guides(fill = guide_legend(nrow = 2))

png('binned_gammadfe_allregions_nes.png', width=6, height=3, units='in', res=320)

dfe_plot

dev.off()
