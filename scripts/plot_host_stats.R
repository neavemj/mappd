
# script to plot the summarised blast files for the subassembly
# to identify the host

library(ggplot2)
library(scales)


plot_host_stats <- function(host_summary, pdf_file, png_file) {
  
  #host_summary <- "/flush3/nea040/mappd_master/sheep/02_depletion/host/tmp.long"
  
  summary_df = read.csv(host_summary, sep="\t", header=T)
  
  # plot the summary table
  p <- ggplot(summary_df, aes(x=reorder(Species, Count), y=Count, fill=Sample)) +
    geom_bar(stat='identity', position='dodge') +
    #theme(axis.text.x = element_text(angle=45, hjust=1)) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    scale_y_continuous(labels = comma) +
    ylab("Reads mapped") +
    #theme(legend.position = "none") +
    coord_flip()
    #facet_wrap(~Sample, scales="free", ncol=1)
  
  # dynamically change figure height depending on number of samples
  # add 2 inches for every additional sample
  ht = 2 * (length(unique(summary_df$Sample)))
  ggsave(pdf_file, p, width=8, height=ht)
  ggsave(png_file, width=8, height=ht, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
host_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_host_stats(host_summary, pdf_file, png_file)


