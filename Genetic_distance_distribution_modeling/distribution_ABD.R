#!/usr/bin/env Rscript

## libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "../log_scale_dist/DD_HH_sample_by_count",
              help="input file"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "DDHH_ABHH_ABDT.pdf",
              help = "[opt] output file name."),
  make_option(c("-n", "--npeak"), dest = "npeak", default = 1,
              help = "[opt] output file name.")
)

parser <- OptionParser(usage = "mapdrawer [options]",
                       option_list = option_list)

## check arguments
arguments <- parse_args(parser)
infile <- arguments$infile
outfile <- arguments$outfile
npeak <- arguments$npeak

#### part1 ----

my_confidence_intervals <- function(DF){
  sample_mean <- mean(DF)
  standard_dev <- sd(DF)
  mean_minus_margin_of_error <- sample_mean - standard_dev
  mean_plus_margin_of_error <- sample_mean + standard_dev
  c(mean_minus_margin_of_error, sample_mean, mean_plus_margin_of_error)
}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma) * 1000
}

turn_label_log <- function(x) {10^-(6-x/10)}
turn_label_log2 <- function(x) {10^-(6-x/10)/1.6*1e8/2}

plot_dist <- function(infile, n, id_thr, lthr, uthr, ratio=1, col1="#3FB1FF", col2="#3F49FF"){
  # n: number of peaks
  # id_thr: ideal threshold, in years
  # lthr: lower threshold, in density *10^6
  # uthr: upper threshold, in density *10^6
  # ratio: scaling factor of lambda
  
  DF1 <- read.table(paste0(infile, ".txt"))
  
  tmp <- data.frame()
  for (i in seq(1,100)){
    t <- my_confidence_intervals(DF1[,i])
    tmp <- rbind(tmp, c(i, t))
  }
  colnames(tmp) <- c("density", "lower", "mean", "upper")
  tmp[tmp<0] <- 0
  
  DF1 <- read.table(paste0(infile, "_flat.txt"))
  
  if (n==1) {
    # library(mclust)
    # mixmdl <- Mclust(DF1$V1, G=1)
    # thr1 <- mixmdl$parameters$mean
    # sig1 <- mixmdl$parameters$variance$sigmasq[1]^0.5
    # lam1 <- mixmdl$parameters$pro[1]
    mixmdl <- mixtools::normalmixEM(DF1$V1, k = 2)
    thr1 <- mixmdl$mu[1]
    print(thr1)
    sig1 <- mixmdl$sigma[1]
    print(sig1)
    lam1 <- mixmdl$lambda[1]
  } else {
    mixmdl <- mixtools::normalmixEM(DF1$V1, k = 2)
    thr1 <- mixmdl$mu[1]
    print(thr1)
    thr2 <- mixmdl$mu[2]
    sig1 <- mixmdl$sigma[1]
    print(sig1)
    sig2 <- mixmdl$sigma[2]
    lam1 <- mixmdl$lambda[1]
    lam2 <- mixmdl$lambda[2]
  }
  
  id_thr <- log10(id_thr*1.6*2/100)*10
  turn_label_log3 <- function(x) {x/1.6/2*100}
  
  lthr <- thr1-sig1
  uthr <- thr1+sig1
  
  # lthr <- log10(lthr)*10
  # uthr <- log10(uthr)*10
  
  p <- ggplot(tmp)+
    geom_ribbon(aes(density, ymin=lower,ymax=upper), fill="grey80", col=NA) +
    annotate(geom = "rect", xmin=lthr, xmax=uthr, ymin=0, ymax=max(tmp$upper), col=NA, fill=alpha("darkgrey", 0.5)) +
    # geom_line(aes(density, mean), size=1.2, col="red") +
    geom_vline(xintercept = 30, col="red", size=2) +
    geom_vline(xintercept = c(thr1), col="black", linetype="longdash") +
    annotate(geom = "text", x=thr1, y=300, label = as.integer(turn_label_log2(thr1)), size=10) +
    annotate(geom = "text", x=15, y=200, label = as.integer(turn_label_log2(lthr)), size=10) +
    annotate(geom = "text", x=15, y=300, label = as.integer(turn_label_log2(uthr)), size=10) +
    stat_function(geom = "line", fun = plot_mix_comps, args = list(thr1, sig1, lam1*ratio), col = col1, lwd = 1.5) +
    # scale_x_continuous(limits=c(0,40), labels = c(expression(10^-6), expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2))) +
    scale_x_continuous(limits=c(10,40), labels = c(expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2))) +
    cowplot::theme_cowplot() +
    theme(text=element_text(family ="Arial")) +
    ylab(NULL) + xlab(NULL)
  
  if (n==2) {
    p <- p + stat_function(geom = "line", fun = plot_mix_comps, args = list(thr2, sig2, lam2*ratio), col = col2, lwd = 1.5) +
      geom_vline(xintercept = c(thr2), col="black", linetype="longdash") +
      annotate(geom = "text", x=thr2, y=300, label = as.integer(turn_label_log2(thr2)), size=10)
  }
  p
}

p1 <- plot_dist("../log_scale_dist/DD_HH_sample_by_count", 1, 6300, 0.7, col1="#FFB111", lthr=184, uthr=215)
p2 <- plot_dist("../log_scale_dist/AABB_HH_sample_by_count", 2, 6300, 2, col1="#3FB1FF", col2="#3F49FF", lthr=151, uthr=262)
p3 <- plot_dist("../log_scale_dist/AABB_DT_sample_by_count", 2, 7100, 2, col1="#3FB1FF", col2="#3F49FF", lthr=153, uthr=334)

p <- cowplot::plot_grid(p1, p2, p3, align = "h", ncol = 1)
ggsave("DDHH_ABHH_ABDT.pdf", p, width = 4, height = 3)

p <- plot_dist(infile, npeak, 6300)
ggsave(outfile, p, width = 5, height = 3)
