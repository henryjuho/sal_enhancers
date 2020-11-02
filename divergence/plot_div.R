library(ggplot2)

setwd('~/sal_enhancers/divergence')
div <- read.delim('all_regs_div.txt')

div_plot <- ggplot(div, aes(x=region, y=divergence)) +
  geom_bar(stat='identity') + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  xlab('')

png('divergence.png', width=3, height=3, res=320, units='in')

div_plot

dev.off()

