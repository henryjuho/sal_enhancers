library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

rm(list=ls())

setwd('~/sal_enhancers/summary_stats')

stats <- read.csv('summary_stats_usfs.csv')

stats$region <- factor(stats$region, levels=rev(c('4fold', 'all_enhancers', 'intergenic', 'intergenic_enhancers', 'intron', 'intron_enhancers', 'utr', 'utr_enhancers', 'cds', 'cds_enhancers', '0fold', '0fold_enhancers')))
levels(stats$region) <- rev(c('4fold degenerate', 'H3K27ac peaks (all)', 'intergenic', 'H3K27ac peaks (intergenic)', 'intron', 'H3K27ac peaks (intron)', 'UTR', 'H3K27ac peaks (UTR)', 'CDS', 'H3K27ac peaks (CDS)', '0fold degenerate', 'H3K27ac peaks (0fold)'))

raw <- subset(stats, bs==0)
bs <- subset(stats, bs!=0)

cis <- summarise(group_by(bs, region),
                 t_lwr=quantile(theta, 0.025), t_upr=quantile(theta, 0.975),
                 p_lwr=quantile(pi, 0.025), p_upr=quantile(pi, 0.975),
                 d_lwr=quantile(tajimas_d, 0.025), d_upr=quantile(tajimas_d, 0.975))

plot_data <- full_join(raw, cis)

str(plot_data)

theta_plot <- ggplot(plot_data, aes(x=region, y=theta)) +
  geom_point(stat='identity', position='dodge') +
  scale_fill_manual(values=viridis(2)) +
  geom_errorbar(aes(ymin=t_lwr, ymax=t_upr), position = position_dodge(width=0.9), width=0.5) +
  theme_classic() +
  labs(x='', y=expression(theta)) +
  theme(legend.title=element_text(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(colour = "black",size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8, angle=45, hjust=1),
        strip.background = element_blank())


pi_plot <- ggplot(plot_data, aes(x=region, y=pi)) +
  geom_point(stat='identity', position='dodge') +
  scale_fill_manual(values=viridis(2)) +
  geom_errorbar(aes(ymin=p_lwr, ymax=p_upr), position = position_dodge(width=0.9), width=0.25) +
  theme_classic() +
  labs(x='', y=expression(pi)) +
  theme(legend.title=element_text(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(colour = "black",size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        strip.background = element_blank()) + coord_flip()

d_plot <- ggplot(plot_data, aes(x=region, y=tajimas_d)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=viridis(2)) +
  geom_errorbar(aes(ymin=d_lwr, ymax=d_upr), position = position_dodge(width=0.9), width=0.25) +
  theme_classic() +
  labs(x='', y="Tajima's D") +
  theme(legend.title=element_text(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        axis.text.y = element_blank(),
        strip.background = element_blank()) + coord_flip()

pi_dat <-subset(plot_data, select=c('region', 'pi', 'p_lwr', 'p_upr'))
colnames(pi_dat) <- c('region', 'value', 'lwr', 'upr')
pi_dat$statistic <- 'pi'

d_dat <-subset(plot_data, select=c('region', 'tajimas_d', 'd_lwr', 'd_upr'))
colnames(d_dat) <- c('region', 'value', 'lwr', 'upr')
d_dat$statistic <- "D"

facet_dat <- rbind(pi_dat, d_dat)

facet_dat$statistic <- factor(facet_dat$statistic, levels=c('pi', "D"), labels=c(expression(pi), expression(Tajimas ~ D)))
vline_dat <- data.frame(statistic=c('pi', "D"), z=c(NA, 0.0))
vline_dat$statistic <- factor(vline_dat$statistic, levels=c('pi', "D"), labels=c(expression(pi), expression(Tajimas ~ D)))

vline_dat

facet_plot <- ggplot(facet_dat, aes(x=region, y=value)) +
  geom_hline(data=vline_dat, aes(yintercept=z), linetype='dashed', colour='grey') +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.25) +
  geom_point(stat='identity', position=position_dodge(width=0.9), shape=21, fill='grey', colour='black') +
  theme_classic() + labs(x='', y='') +
  facet_wrap(~statistic, nrow=1, scales='free_x', strip.position="bottom", labeller = label_parsed) +
  theme(legend.title=element_text(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(colour = "black",size=8),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_text(colour = "black",size=8),
        strip.background = element_blank(), strip.placement = "outside") + coord_flip()



layout <- rbind(c(1,1,1,1,2,2,2))

png('summary_stats.png', res=320, units='in', height=2.5, width=6)

facet_plot
#grid.arrange(pi_plot, d_plot, layout_matrix = layout)

dev.off()

pdf('summary_stats.pdf', height=2.5, width=6)

facet_plot
#grid.arrange(pi_plot, d_plot, layout_matrix = layout)

dev.off()