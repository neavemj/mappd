
# script to plot the summarised log files from trimmomatic

library(ggplot2)
library(tidyr)
require(scales)


plot_trim <- function(trim_summary, pdf_file, png_file) {
  
  summary_df = read.table(trim_summary, sep="\t", header=T)

  summary_long = gather(summary_df, Result, pairs, both_surviving, forward_only, reverse_only, dropped)

  summary_long$Result <- factor(summary_long$Result, levels = c("dropped", "reverse_only", "forward_only", "both_surviving"))

  cols <- c("both_surviving" = "#7CAE00", "forward_only" = "#C77CFF", "reverse_only" = "#00BFC4", "dropped" = "#F8766D")
  
  p <- ggplot(summary_long, aes(x=sample, y=pairs, fill=Result)) +
    geom_bar(stat='identity') +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values=cols) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ylab("Read pairs") +
    xlab("Sample")
  
  # make the figure width adjust dynamically with the number of samples
  # add 1/3 of an inch for each additional samples
  # TODO: test with say 10 or 20 samples
  fig_width = 4 + (length(unique(summary_long$sample)) / 3)
   
  ggsave(pdf_file, p, width=fig_width, height=4)
  #ggsave(png_file, width=fig_width, height=4, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
trim_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_trim(trim_summary, pdf_file, png_file)


