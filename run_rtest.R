args = commandArgs(trailingOnly=TRUE)
library(CC2Sim)

m <- args[1]
dat <- read.csv(args[2])
out_path <- args[3]

methods = list(
  DESeq2 = run.DESeq2,
  edgeR = run.edgeR,
  sgRSEA = run.sgRSEA,
  PBNPA = run.PBNPA
)

ret <- methods[[m]](dat)

sgrna_path <- file.path(out_path, "sgRNA.csv")
gene_path <- file.path(out_path, "gene.csv")

if(!is.null(ret$sgRNA)) write.csv(ret$sgRNA, sgrna_path)
if(!is.null(ret$gene)) write.csv(ret$gene, gene_path)
