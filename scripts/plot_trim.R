
# script to plot the summarised log files from trimmomatic

library(ggplot2)
library(tidyr)
require(scales)


plot_trim <- function(trim_summary, pdf_file, png_file) {
  
  #if (makePDF) {
    #CREATE pdf as output file
  #  pdf(file = pdf_file)
  #} else {
    #CREATE png as output file
  #  png(file=png_file, width = 8, height = 8, unit="in",res=300)
  #}
  
  summary_df = read.table(trim_summary, sep="\t", header=T)

  summary_long = gather(summary_df, type, pairs, both_surviving, forward_only, reverse_only, dropped)

  summary_long$type <- factor(summary_long$type, levels = c("dropped", "reverse_only", "forward_only", "both_surviving"))

  cols <- c("both_surviving" = "#7CAE00", "forward_only" = "#C77CFF", "reverse_only" = "#00BFC4", "dropped" = "#F8766D")
  
  ggplot(summary_long, aes(x=sample, y=pairs, fill=type)) +
    geom_bar(stat='identity') +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values=cols) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  ggsave(pdf_file, width=4, height=4)
  ggsave(png_file, width=4, height=4, dpi=300)

  
}

args <- commandArgs(trailingOnly = TRUE)
trim_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_trim(trim_summary, pdf_file, png_file)


