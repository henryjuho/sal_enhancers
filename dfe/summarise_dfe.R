library(ggplot2)
library(dplyr)
library(viridis)
library(gridExtra)
library(comprehenr)

setwd('~/sal_enhancers/dfe')

c1 <- viridis(n=6)
c2 <- adjustcolor(c1, alpha.f=0.5)
colours <- to_vec(for(i in 1:6) c(c1[i], c2[i]))

#=======================================================================================================================
#     DFE
#=======================================================================================================================

# dfe <- read.csv('binned_dfe_allregions.csv')
dfe <- read.csv('binned_gammadfe_allregions_nes.csv')

str(dfe)
dfe$region <- factor(dfe$region, levels=c('intergenic', 'intergenic-enhancers', 'intron', 'intron-enhancers', 'utr', 'utr-enhancers', 'cds', 'cds-enhancers', '0fold', '0fold-enhancers', 'all-enhancers'))
levels(dfe$region) <- c('intergenic', 'H3K27ac peaks (intergenic)', 'intron', 'H3K27ac peaks (intron)', 'UTR', 'H3K27ac peaks (UTR)', 'CDS', 'H3K27ac peaks (CDS)', '0fold', 'H3K27ac peaks (0fold)', 'H3K27ac peaks (all)')

levels(dfe$bin) <- c("0 - 1\nneutral", "1 - 10\nweakly deleterious", "10 - 100\ndeleterious", "> 100\nstrongly deleterious")
str(dfe)

raw <- subset(dfe, rep==0)
bs <- subset(dfe, rep!=0)

cis <- summarise(group_by(bs, region, bin),
                 lwr=quantile(proportion, 0.025), upr=quantile(proportion, 0.975))

plot_data <- full_join(raw, cis)

dfe_plot <- ggplot(plot_data, aes(x=bin, y=proportion, fill=region)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_manual(values=colours) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.5) +
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
  xlab(expression(N[e] * s)) + ylab('proportion of variants') +
  ggtitle('(a)') +
  guides(fill = guide_legend(nrow = 2))

#=======================================================================================================================
#     Alpha
#=======================================================================================================================

alp <- read.csv('../divergence/alpha_estimates.csv')

str(alp)
alp$region <- factor(alp$region, levels=c('intergenic', 'intergenic-enhancers', 'intron', 'intron-enhancers', 'utr', 'utr-enhancers', 'cds', 'cds-enhancers', '0fold', '0fold-enhancers', 'all-enhancers'))
levels(alp$region) <- c('intergenic', 'H3K27ac peaks (intergenic)', 'intron', 'H3K27ac peaks (intron)', 'UTR', 'H3K27ac peaks (UTR)', 'CDS', 'H3K27ac peaks (CDS)', '0fold', 'H3K27ac peaks (0fold)', 'H3K27ac peaks (all)')

alp$bin <- 'alpha'

raw <- subset(alp, bs_rep==0)
bs <- subset(alp, bs_rep!=0)

cis <- summarise(group_by(bs, region),
                 lwr=quantile(alpha, 0.025), upr=quantile(alpha, 0.975))

plot_data <- full_join(raw, cis)

alpha_plot <- ggplot(plot_data, aes(x=bin, y=alpha, fill=region)) +
  geom_bar(stat='identity', position = position_dodge(width=0.9)) +
  scale_fill_manual(values=colours) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position = position_dodge(width=0.9), width=0.5) +
  theme_classic() + labs(x='', y=expression(alpha)) +
  theme(legend.title=element_blank(),
        legend.position = 'none',
        legend.margin = margin(t = -1, r = 0, b = 2, l = -10, unit = "mm"),
        legend.text = element_text(colour = "black",size=6),
        axis.text.y = element_text(colour = "black",size=8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=9, vjust=-6, hjust=0.02, face='bold'),
        plot.margin = margin(t=0, r=5.5, b=24, l=0, unit = "pt"))+
  ggtitle('(b)') +
  guides(fill = guide_legend(nrow = 2))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


leg <- g_legend(dfe_plot)

png('binned_gammadfe_allregions_nes.png', width=7.5, height=3, units='in', res=320)

layout <- rbind(c(1,1,1,1,1,1,2,2),
                c(3,3,3,3,3,3,3,3))

grid.arrange(dfe_plot + theme(legend.position="none"), alpha_plot, leg, layout_matrix=layout, heights=c(10,1))

dev.off()

pdf('binned_gammadfe_allregions_nes.pdf', width=7.5, height=3)

layout <- rbind(c(1,1,1,1,1,1,2,2),
                c(3,3,3,3,3,3,3,3))

grid.arrange(dfe_plot + theme(legend.position="none"), alpha_plot, leg, layout_matrix=layout, heights=c(10,1))

dev.off()