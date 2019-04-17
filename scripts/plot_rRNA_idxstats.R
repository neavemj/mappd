
# script to plot the summarised count files for the rRNA databases

library(ggplot2)
library(tidyverse)
library(scales)


plot_idxstats <- function(idxstats_summary, pdf_file, png_file, tsv_file) {
  
  #idxstats_summary <- "/flush3/nea040/mappd_master/freshwater_prawn_mappd/02_depletion/rRNA/idxstats.summary"
  
  summary_df = read.csv(idxstats_summary, sep="\t", header=T)
  
  # decided not to plot or summarise the LSU results 
  summary_df = subset(summary_df, rRNA_Type != "LSU")
  
  species_summary <- summary_df %>%
    group_by(Sample, Species) %>%
    summarise(Sum_Mapped = sum(Mapped_Reads)) %>%
    top_n(10, Sum_Mapped)

  # write out a summary table for the report
  # need to add the taxonomy string to the summary file
  # this returns just the first match for taxonomy string for each species (called 'Organism')
  species_summary_tax <- merge(species_summary, aggregate(Species ~ Organism, data=summary_df, head, 1), by="Species") 
  
  # also add taxid
  species_summary_taxid <- merge(species_summary_tax, aggregate(Species ~ Taxid, data=summary_df, head, 1), by="Species")
  
  # change the data to wide format for a nicer table
  species_summary_tax_long <- spread(species_summary_taxid[,c("Sample", "Sum_Mapped", "Organism", "Taxid")], Sample, Sum_Mapped, fill=0)
  
  # sort by the most abundant species (first sample column only)
  species_summary_tax_long <- species_summary_tax_long[order(-species_summary_tax_long[,3]), ]
  
  # write summary table to include in the report later
  write.table(species_summary_tax_long, tsv_file, quote=F, sep="\t", row.names = F)
  
  # plot the summary table
  p <- ggplot(species_summary, aes(x=reorder(Species, Sum_Mapped), y=Sum_Mapped, fill=Sample)) +
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
idxstats_summary = args[1]
pdf_file = args[2]
png_file = args[3]
tsv_file = args[4]
plot_idxstats(idxstats_summary, pdf_file, png_file, tsv_file)


