#!/usr/bin/env Rscript
# module load R/R-3.6.3
library("ggplot2")
df <- read.delim(file = "benchmarks2.tsv", sep = '\t')
colnames(df) <- c("num_threads", "time", "num_pairs", "num_tumors", "num_normals", "action")
df[["time_per_pair"]] <- df[["time"]] / df[["num_pairs"]]
df[["num_threads"]] <- factor(df[["num_threads"]], levels = sort(unique(df[["num_threads"]])))

pdf("benchmark_time_pairs_threads.pdf", width = 10)
ggplot(data = df, aes(x = num_pairs, y = time/60, color = num_threads, fill = action)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme_bw() +
    ggtitle("Conpair Execution Time (min) vs. Number of Pairs vs. Thread Count")
dev.off()

pdf("time_per_pair.pdf")
ggplot(data = df, aes(x = num_threads, y = time_per_pair)) +
    geom_boxplot() +
    ggtitle("Conpair Execution Time Per Pair (sec) vs. Thread Count") +
    theme_bw()
dev.off()

save.image()
