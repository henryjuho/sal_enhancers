library(ggplot2)
library(dplyr)
library(viridis)
library(gridExtra)

setwd('~/sal_enhancers/divergence')

alp <- read.csv('alpha_estimates.csv')

str(alp)
alp$region <- factor(alp$region, levels=c('intergenic', 'intergenic-enhancers', 'intron', 'intron-enhancers', 'utr', 'utr-enhancers', 'cds', 'cds-enhancers', 'all-enhancers'))
levels(alp$region) <- c('intergenic', 'H3K27ac peaks (intergenic)', 'intron', 'H3K27ac peaks (intron)', 'UTR', 'H3K27ac peaks (UTR)', 'CDS', 'H3K27ac peaks (CDS)', 'H3K27ac peaks (all)')

alp$bin <- 'alpha'

raw <- subset(alp, bs_rep==0)
bs <- subset(alp, bs_rep!=0)

cis <- summarise(group_by(bs, region),
                 lwr=quantile(alpha, 0.025), upr=quantile(alpha, 0.975))

plot_data <- full_join(raw, cis)

alpha_plot <- ggplot(plot_data, aes(x=bin, y=alpha, fill=region)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_manual(values=viridis(10)[2:10]) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.5) +
  theme_classic() + labs(y='', x=expression(alpha)) +
  theme(legend.title=element_blank(),
        legend.position = 'none',
        legend.margin = margin(t = -1, r = 0, b = 2, l = -10, unit = "mm"),
        legend.text = element_text(colour = "black",size=6),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))

omega_cis <- summarise(group_by(bs, region),
                 lwr=quantile(omega, 0.025), upr=quantile(omega, 0.975))

plot_data <- full_join(raw, omega_cis)

omega_plot <- ggplot(plot_data, aes(x=bin, y=omega, fill=region)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_manual(values=viridis(10)[2:10]) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.5) +
  theme_classic() + labs(y='', x=expression(omega[alpha])) +
  theme(legend.title=element_blank(),
        legend.position = 'none',
        legend.margin = margin(t = -1, r = 0, b = 2, l = -10, unit = "mm"),
        legend.text = element_text(colour = "black",size=6),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_blank(),
        strip.background = element_blank())+
  guides(fill = guide_legend(nrow = 2))


png('alpha.png', width=3, height=3, res=320, units='in')

#grid.arrange(alpha_plot, omega_plot, nrow=1)
alpha_plot

dev.off()