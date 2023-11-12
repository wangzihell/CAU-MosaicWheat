#!/usr/bin/env Rscript

## libraries
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "chr1D.centromo",
              help="input file"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "chr1D.centromo.cluster",
              help = "[opt] output file name."),
  make_option(c("-r", "--ratio"), dest = "ratio", default = 0.5,
              help = "[opt] output file name.") 
)

parser <- OptionParser(usage = "mapdrawer [options]",
                       option_list = option_list)

## check arguments
arguments <- parse_args(parser)
infile <- arguments$infile
outfile <- arguments$outfile
ratio <- as.numeric(arguments$ratio)

DF1 <- read.table(infile, header = T)
samples <- names(DF1)[-c(1,2,3)]
DF2 <- abs(as.matrix(DF1[,-c(1,2,3)]))

n <- ncol(DF2)
m <- nrow(DF2)
dis_matrix <- matrix(nrow = n, ncol = n)
Sum <- 0
for (i in 1:n) {
  for (j in 1:n) {
    Sum <- sum(DF2[,i] != DF2[,j], na.rm = T)
    dis_matrix[i,j] <- Sum
    Sum <- 0
  }
}
rownames(dis_matrix) <- samples
colnames(dis_matrix) <- samples
DF3 <- cutree(hclust(as.dist(dis_matrix)), h = m*ratio)
write.table(DF3, outfile, row.names = T, col.names = F, quote = F, sep = "\t")
