
# script to plot the summarised count files for the rRNA databases

library(ggplot2)
library(dplyr)
library(scales)


plot_idxstats <- function(idxstats_summary, pdf_file, png_file) {
  
  #idxstats_summary <- "/flush3/nea040/mappd_master/freshwater_prawn_mappd/02_depletion/LSU.idxstats.summary"
  
  summary_df = read.csv(idxstats_summary, sep="\t", header=T)
  
  species_summary <- summary_df %>%
    group_by(Sample, Species) %>%
    summarise(Sum_Mapped = sum(Mapped_Reads)) %>%
    top_n(10, Sum_Mapped)


  p <- ggplot(species_summary, aes(x=reorder(Species, Sum_Mapped), y=Sum_Mapped)) +
    geom_bar(stat='identity') +
    #theme(axis.text.x = element_text(angle=45, hjust=1)) +
    theme(axis.title.y = element_blank()) +
    scale_y_continuous(labels = comma) +
    ylab("Reads mapped") +
    theme(legend.position = "none") + 
    coord_flip() +
    facet_wrap(~Sample, scales="free", ncol=1)
  
  # dynamically change figure height depending on number of samples
  ht = length(unique(summary_df$Sample)) * 3
  ggsave(pdf_file, p, width=6, height=ht)
  ggsave(png_file, width=6, height=ht, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
idxstats_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_idxstats(idxstats_summary, pdf_file, png_file)


