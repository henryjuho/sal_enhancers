library(ggplot2)
library(dplyr)
library(viridis)


dfe <- read.csv('binned_dfe_allregions.csv')
dfe$region <- factor(dfe$region, levels=c('intron', 'utr', 'cds'))

raw <- subset(dfe, rep==0)
bs <- subset(dfe, rep!=0)

cis <- summarise(group_by(bs, region, bin),
                 lwr=quantile(proportion, 0.025), upr=quantile(proportion, 0.975))

plot_data <- full_join(raw, cis)

dfe_plot <- ggplot(plot_data, aes(x=bin, y=proportion, fill=region)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_manual(values=viridis(3)) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.5) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.position=c(0.8, 0.75)) +
  xlab(expression(gamma == 4 * N[e] * s))

png('all_regions_dfe.png', width=3, height=2.5, units='in', res=320)

dfe_plot

dev.off()
