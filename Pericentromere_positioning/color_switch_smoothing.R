#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rlang))
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = NULL, 
              help = "the chromosome name, such as chr1A"),
  make_option(c("-c", "--callback"), dest = "callback", default = "mean", 
              help = "the chromosome name, such as chr1A"),
  make_option(c("-n", "--n_col"), dest = "n_col", default = 4, 
              help = "the chromosome name, such as chr1A"),
  make_option(c("-w", "--width"), dest = "width", default = 4,
              help = "the chromosome name, such as chr1A"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = NULL, 
              help = "the chromosome name, such as chr1A")
)
#
parser <- OptionParser(usage = "%prog [options] file", 
                       option_list=option_list)

#
arguments <- parse_args(parser)


DF1 <- read.table(arguments$infile)
n_col <- arguments$n_col
width <- arguments$width
outfile <- arguments$outfile
callback <- arguments$callback

sliding_window <- function(x,              # numeric vector
                           width = 8,          # window size
                           callback = "mean",  # summary function to return a single value
                           fill = F) {   # make output as long as input? Set TRUE for adding data.frame columns
  
  x_length <- length(x)
  stopifnot(x_length > width)  # check if this is appropriate
  
  starts <- seq_len(x_length - width)  # start points
  ends <- starts + width  # end points
  
  if ( fill ) {  # make adjustments if output length should match input length
    starts <- seq_len(x_length)
    ends <- c(ends, rep(x_length, width))
  }
  
  tuples <- mapply(c, starts, ends, SIMPLIFY=FALSE)
  if (callback=="mean"){
    vector_f <- function(y) mean(x[y[1]:y[2]])
  } else {
    vector_f <- function(y) median(x[y[1]:y[2]])
  }
  results <- vapply(tuples, vector_f, vector_f(tuples[[1]]))  # vapply gives a speed-up (in theory)
  return(c(rep(median(x[1:width/2]), width/2), results, rep(median(x[x_length - width/2:x_length]), width/2)))
}

DF1[,n_col] <- sliding_window(DF1[,n_col], width=width, callback=callback)
write.table(DF1, outfile, row.names = F, col.names = F, quote = F, sep = "\t")