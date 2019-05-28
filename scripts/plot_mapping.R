
# script to plot the summarised mapping files from bowtie

library(ggplot2)
library(tidyr)
require(scales)


plot_mapping <- function(mapping_summary, pdf_file, png_file) {

  #mapping_summary <- "/flush3/nea040/mappd_master/freshwater_prawn_mappd/logs/rRNA_mapping_summary.tsv"

  summary_df = read.table(mapping_summary, sep="\t", header=T)

  summary_df$Type <- factor(summary_df$Type, levels = c("rRNA_LSU", "rRNA_SSU", "host", "surviving"))

  cols <- c("surviving" = "#7CAE00", "rRNA_LSU" = "#C77CFF", "rRNA_SSU" = "#00BFC4", "host" = "#F8766D")

  p <- ggplot(summary_df, aes(x=Sample, y=Paired_Reads, fill=Type)) +
    geom_bar(stat='identity') +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values=cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ylab("Read pairs") +
    xlab("Sample")


  # make the figure width adjust dynamically with the number of samples
  # add 1/3 of an inch for each additional samples
  # TODO: test with say 10 or 20 samples
  fig_width = 4 + (length(unique(summary_df$Sample)) / 3)

  ggsave(pdf_file, p, width=fig_width, height=4)
  ggsave(png_file, width=fig_width, height=4, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
mapping_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_mapping(mapping_summary, pdf_file, png_file)


