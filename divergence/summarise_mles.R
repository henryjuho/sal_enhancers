library(dplyr)

setwd('~/sal_enhancers/divergence')

mles <- read.csv('alpha_estimates.csv')

mles$region <- factor(mles$region, levels=rev(c('all-enhancers', 'intergenic-enhancers', 'intergenic', 'intron-enhancers', 'intron','utr-enhancers',  'utr', 'cds-enhancers', 'cds', '0fold')))
levels(mles$region) <- rev(c( 'H3K27ac peaks (all)', 'H3K27ac peaks (intergenic)', 'intergenic', 'H3K27ac peaks (intron)', 'intron', 'H3K27ac peaks (UTR)', 'UTR', 'H3K27ac peaks (CDS)', 'CDS', '0fold degenerate'))

raw <- subset(mles, bs_rep==0)
bs <- subset(mles, bs_rep!=0)

cis <- summarise(group_by(bs, region),
                 t_lwr=quantile(theta, 0.025), t_upr=quantile(theta, 0.975),
                 sh_lwr=quantile(shape, 0.025), sh_upr=quantile(shape, 0.975),
                 sc_lwr=quantile(scale, 0.025), sc_upr=quantile(scale, 0.975),
                 e_lwr=quantile(error, 0.025), e_upr=quantile(error, 0.975),
                 a_lwr=quantile(alpha, 0.025), a_upr=quantile(alpha, 0.975))

sum_dat <- dplyr::full_join(raw, cis)

out_dat <- data.frame(region=sum_dat$region,
                      theta=paste(signif(sum_dat$theta, 3), ' (', signif(sum_dat$t_lwr, 3), ' ', signif(sum_dat$t_upr, 3), ')', sep=''),
                      shape=paste(signif(sum_dat$shape, 3), ' (', signif(sum_dat$sh_lwr, 3), ' ', signif(sum_dat$sh_upr, 3), ')', sep=''),
                      scale=paste(signif(sum_dat$scale, 3), ' (', signif(sum_dat$sc_lwr, 3), ' ', signif(sum_dat$sc_upr, 3), ')', sep=''),
                      error=paste(signif(sum_dat$error, 3), ' (', signif(sum_dat$e_lwr, 3), ' ', signif(sum_dat$e_upr, 3), ')', sep=''),
                      alpha=paste(signif(sum_dat$alpha*100, 3), ' (', signif(sum_dat$a_lwr*100, 3), ' ', signif(sum_dat$a_upr*100, 3), ')', sep=''))

out_dat <- with(out_dat, out_dat[order(region),])

write.csv(out_dat, 'mle_summary.csv', quote=FALSE, row.names=FALSE)