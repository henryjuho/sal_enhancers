library(ape)

args = commandArgs(TRUE)

fasta = args[1]
region = args[2]

data = read.FASTA(fasta)

divergance = dist.dna(data, model = "K80")

sal_salb = divergance[1]
sal_pike = divergance[2]
salb_pike = divergance[3]

sal = (sal_pike - salb_pike + sal_salb) / 2
salb = (salb_pike - sal_pike + sal_salb) / 2

cat(region, sal_salb, sal, salb, sep=',')
