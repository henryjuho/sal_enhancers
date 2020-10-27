library(ape)

args <- commandArgs(TRUE)

fasta <- args[1]
region <- args[2]

data <- read.FASTA(fasta)

divergance <- dist.dna(data, model = "K80")

cat(paste('region', 'divergence', sep='\t'))
cat('\n')

sal_charr <- divergance[1]
sal_trout <- divergance[2]
charr_trout <- divergance[3]

sal <- (sal_charr - charr_trout + sal_trout) / 2

cat(region, sal, sep='\t')
cat('\n')
