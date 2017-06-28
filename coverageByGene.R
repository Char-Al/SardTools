#!/usr/bin/env Rscript

library(reshape2)
library(ggplot2)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "./"
}

# raw <- read.csv("sample_DPamplicons.csv",sep=";")
raw <- read.csv(args[1],sep=";")
raw$gene <- gsub("Ex.+$", "", raw$PCR_name)
raw.m <- melt(raw,id.vars=c('sample','PCR_name','gene'))

g <- ggplot(raw.m, aes(x = PCR_name, y = value, color = variable, group = variable))
g <- g + geom_point() +
	geom_line() +
	scale_y_log10() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	facet_wrap(~gene, scales = "free_x")

ggsave(filename=paste0(args[2],"/plot.pdf"), plot=g, width = 21, height = 29, units = "cm")
