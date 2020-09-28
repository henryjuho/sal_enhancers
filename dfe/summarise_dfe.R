library(ggplot2)
library(dplyr)
library(viridis)

setwd('/Users/henryjuho/sal_enhancers/dfe')

dfe <- read.csv('binned_dfe_allregions.csv')
str(dfe)
dfe$region <- factor(dfe$region, levels=c('intergenic-enhancers', 'intron', 'intron-enhancers', 'utr', 'utr-enhancers', 'cds'))

str(dfe)

raw <- subset(dfe, rep==0)
bs <- subset(dfe, rep!=0)

cis <- summarise(group_by(bs, region, bin),
                 lwr=quantile(proportion, 0.025), upr=quantile(proportion, 0.975))

plot_data <- full_join(raw, cis)

dfe_plot <- ggplot(plot_data, aes(x=bin, y=proportion, fill=region)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_manual(values=viridis(7)[2:7]) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.5) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = "bottom", legend.key.size = unit(.33, "cm"),
        legend.text = element_text(colour = "black",size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        strip.background = element_blank())+
  xlab(expression(gamma == 4 * N[e] * s)) +
  guides(fill = guide_legend(nrow = 2))

png('all_regions_dfe.png', width=4, height=3, units='in', res=320)

dfe_plot

dev.off()
