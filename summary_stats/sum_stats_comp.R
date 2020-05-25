library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

rm(list=ls())

setwd('/Users/henryjuho/sal_enhancers/summary_stats')

sel_dat <- read.csv('summary_stats_sel_sites_usfs.csv')
neu_dat <- read.csv('summary_stats_neu_sites_usfs.csv')

stats <- rbind(sel_dat, neu_dat)

stats$region <- factor(stats$region, levels=c('4fold', 'intergenic_enhancers', 'intron_enhancers', 'utr_enhancers', 'intron', 'utr', 'cds'))

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
  theme_bw() +
  labs(x='', y=expression(theta)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

pi_plot <- ggplot(plot_data, aes(x=region, y=pi)) +
  geom_point(stat='identity', position='dodge') +
  scale_fill_manual(values=viridis(2)) +
  geom_errorbar(aes(ymin=p_lwr, ymax=p_upr), position = position_dodge(width=0.9), width=0.5) +
  theme_bw() +
  labs(x='', y=expression(pi)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

d_plot <- ggplot(plot_data, aes(x=region, y=tajimas_d)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=viridis(2)) +
  geom_errorbar(aes(ymin=d_lwr, ymax=d_upr), position = position_dodge(width=0.9), width=0.5) +
  theme_bw() +
  labs(x='', y="Tajima's D") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

png('summary_stats.png', res=320, units='in', height=3, width=9)

grid.arrange(theta_plot, pi_plot, d_plot, nrow=1)

dev.off()
