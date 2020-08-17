#!/usr/bin/env Rscript

# script to plot some concordance metrics

# module load R/R-3.6.3
library("ggplot2")

args <- commandArgs(TRUE)
input_file <- args[1]
df <- read.delim(file = input_file, header = TRUE, sep = '\t', stringsAsFactors = TRUE)

# drop some columns for faster processing..
# keep_cols <- c("concordance", "tumor", "normal")
# df <- df[, colnames(df)[colnames(df) %in% keep_cols]]

# sort the tumor factor levels for plotting
df[["tumor"]] <- factor(x = df[["tumor"]], levels = sort(levels(df[["tumor"]])))
df[["normal"]] <- factor(x = df[["normal"]], levels = sort(levels(df[["normal"]])))
df[["pcnt_markers"]] <- df[["num_markers_used"]] / df[["num_total_markers"]]
df[["confidence"]] <- df[["concordance"]] * df[["pcnt_markers"]]
save.image()

# # default height
fig_height <- 12
# recalculate for >20 samples
num_tumors <- length(levels(df[["tumor"]]))
if ( num_tumors > 20) {
    fig_height <-  num_tumors / 2
}


pdf(file = "concordance_dist.pdf", height = fig_height)
ggplot(data = df, aes(y = concordance, x = tumor)) +
    geom_violin() +
    coord_cartesian(ylim = c(0, 1)) + # make sure y scale goes the full distribution
    coord_flip() +
    theme_bw() +
    ggtitle("Distribution of concordances per tumor")
dev.off()

pdf(file = "pcnt_markers_dist.pdf", height = fig_height)
ggplot(data = df, aes(y = pcnt_markers, x = tumor)) +
    geom_violin() +
    coord_cartesian(ylim = c(0, 1)) + # make sure y scale goes the full distribution
    coord_flip() +
    theme_bw() +
    ggtitle("Distribution of percent of markers used per tumor")
dev.off()


# df[["normal"]]
# summary(df)
# plot(df[["concordance"]] / df[["num_markers_used"]])
# box
