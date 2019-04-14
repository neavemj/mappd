
# script to plot the summarised count files for the rRNA databases

library(ggplot2)
library(dplyr)
library(scales)


plot_idxstats <- function(idxstats_summary, pdf_file, png_file) {
  
  #idxstats_summary <- "/flush3/nea040/mappd_master/freshwater_prawn_mappd/02_depletion/prawn1_LSU.idxstats.summary"
  
  summary_df = read.table(idxstats_summary, sep="\t", header=T)
  
  species_summary <- summary_df %>%
    group_by(species) %>%
    summarise(sum_mapped = sum(mapped_reads)) %>%
    top_n(10, sum_mapped)


  p <- ggplot(species_summary, aes(x=reorder(species, -sum_mapped), y=sum_mapped)) +
    geom_bar(stat='identity') +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_y_continuous(labels = comma) +
    ylab("Reads mapped") +
    xlab("Species") +
    theme(legend.position = "none")
   
  ggsave(pdf_file, p, width=8, height=5)
  ggsave(png_file, width=8, height=5, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
idxstats_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_idxstats(idxstats_summary, pdf_file, png_file)


