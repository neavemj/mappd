
# script to plot the abundance files for superkingdo

library(ggplot2)
library(tidyverse)
library(scales)
library(RColorBrewer)


plot_abund <- function(abund_summary, tsv_file, pdf_file, png_file) {

  #abund_summary <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/04_annotation/diamond/diamond_blastx_abundance.euk"

  summary_df = read.csv(abund_summary, sep="\t", header=T)

  # the summarise bit keeps family as a column - gets dropped otherwise
  species_summary <- summary_df %>%
    group_by(Sample, Species) %>%
    summarise(Sum_Mapped = sum(Reads_Mapped), Family=toString(Family)) %>%
    top_n(10, wt = Sum_Mapped)
  
  species_summary <- subset(species_summary, Species != "<not present>")
  
  # if a species is not in a particular sample, ggplot draws a wide column
  # want to have a 0 in this case so the column width is constant
  # also want to create a table in wide format
  species_summary_wide = spread(species_summary, Sample, Sum_Mapped, fill=0)
  
  # this calculates rowsums for the samples
  # sorts them by the total and gets just the top 10
  summary_table <- species_summary_wide %>%
    mutate(Total = rowSums(.[3:ncol(species_summary_wide)])) %>%
    arrange(wt = -Total) %>%
    top_n(10, wt = Total)
  
  # write wide format table
  write.table(summary_table, tsv_file, quote=F, sep="\t", row.names = F)
  
  species_summary_long =  gather(summary_table, Sample, Sum_Mapped, 3:ncol(species_summary_wide))
  
  # flip everything around due to my coord_flip() in ggplot call
  species_summary_long$Sample <- as.factor(species_summary_long$Sample)
  species_summary_long$Sample <- factor(species_summary_long$Sample, levels=rev(levels(species_summary_long$Sample)))
  
  # brewer pallete colours - Paired or Set3 for 10 samples max
  cols <- brewer.pal(length(unique(summary_table$Species)), "Set3")
  
  p <- ggplot(species_summary_long, aes(x=Sample, y=Sum_Mapped, fill=Species)) +
    geom_bar(stat='identity') +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values=cols) +
    ylab("Reads") +
    coord_flip()
  
  # dynamically change figure height depending on number of samples
  # although, has to be at least 3 inches high for the legend
  # add 1 inch for every additional sample
  
  if(length(unique(species_summary_long$Sample)) < 3){
    ht = 3
  } else {
    ht = 1 * (length(unique(species_summary_long$Sample)))
  }
  
  ggsave(pdf_file, p, width=8, height=ht)
  ggsave(png_file, width=8, height=ht, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
abund_summary = args[1]
tsv_file = args[2]
pdf_file = args[3]
png_file = args[4]

plot_abund(abund_summary, tsv_file, pdf_file, png_file)


