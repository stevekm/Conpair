#!/usr/bin/env Rscript
# module load R/R-3.6.3
library("ggplot2")
df <- read.delim(file = "benchmarks.tsv", sep = '\t')
colnames(df) <- c("num_threads", "time", "num_pairs", "num_tumors", "num_normals", "action")
df[["time_per_pair"]] <- df[["time"]] / df[["num_pairs"]]
df[["num_threads"]] <- factor(df[["num_threads"]], levels = sort(unique(df[["num_threads"]])))

pdf("time_pairs_threads.pdf", width = 10)
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

aggr <- aggregate(time_per_pair~num_threads, data = df, mean)
colnames(aggr) <- c("num_threads", "mean")
aggr[["median"]] <- aggregate(time_per_pair~num_threads, data = df, median)[["time_per_pair"]]
aggr[["min"]] <- aggregate(time_per_pair~num_threads, data = df, min)[["time_per_pair"]]
aggr[["max"]] <- aggregate(time_per_pair~num_threads, data = df, max)[["time_per_pair"]]
aggr[["sd"]] <- aggregate(time_per_pair~num_threads, data = df, sd)[["time_per_pair"]]

write.table(x = aggr, file = "aggregate_time_per_pair.tsv", quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

save.image()
