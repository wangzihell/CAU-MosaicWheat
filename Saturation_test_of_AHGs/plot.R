#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "saturation_100_LR.percentile.txt",
              help = "dist matrix file, n rows by n cols. [default: whole_all.mdist]"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "4groups_saturation.pdf",
              help = "output prefix.")
)

parser <- OptionParser(usage = "Rscript treeWholeNJ.R [options] Tip group files",
                       option_list = option_list)

arguments <- parse_args(parser)

infile <- arguments$infile
outfile <- arguments$outfile

DF1 <- read.table("saturation_100_WE.percentile.txt")
DF2 <- read.table("saturation_100_DO.percentile.txt")
DF3 <- read.table("saturation_100_LR.percentile.txt")
DF4 <- read.table("saturation_100_CV.percentile.txt")

ggplot() +
  geom_ribbon(data=DF1, aes(x=V1,ymin=V2,ymax=V4), fill=alpha("red", 0.1)) + 
  geom_line(data=DF1, aes(x=V1,y=V3), col="red") +
  geom_ribbon(data=DF2, aes(x=V1,ymin=V2,ymax=V4), fill=alpha("blue", 0.1)) + 
  geom_line(data=DF2, aes(x=V1,y=V3), col="blue") +
  geom_ribbon(data=DF3, aes(x=V1,ymin=V2,ymax=V4), fill=alpha("green", 0.1)) + 
  geom_line(data=DF3, aes(x=V1,y=V3), col="green") +
  geom_ribbon(data=DF4, aes(x=V1,ymin=V2,ymax=V4), fill=alpha("purple", 0.1)) + 
  geom_line(data=DF4, aes(x=V1,y=V3), col="purple") +
  labs(x="Accession number", y="Type number") +
  geom_abline(intercept=c(10150, 16240), color="blue") +
  cowplot::theme_cowplot()
  #geom_point(aes(V1,V2))
ggsave(outfile, width = 2.2, height = 2)


DF1 <- read.table("saturation_100_Western.percentile.txt")
DF2 <- read.table("saturation_100_Eastern.percentile.txt")

ggplot() +
  geom_ribbon(data=DF1, aes(x=V1,ymin=V2,ymax=V4), fill=alpha("red", 0.1)) + 
  geom_line(data=DF1, aes(x=V1,y=V3), col="red") +
  geom_ribbon(data=DF2, aes(x=V1,ymin=V2,ymax=V4), fill=alpha("blue", 0.1)) + 
  geom_line(data=DF2, aes(x=V1,y=V3), col="blue") +
  labs(x="Accession number", y="Type number") +
  cowplot::theme_cowplot()
ggsave("LR_WestEast_saturation.pdf", width = 3, height = 3)
