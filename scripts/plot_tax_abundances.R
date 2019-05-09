
# script to plot the abundance files for superkingdo

library(ggplot2)
library(tidyverse)
library(scales)


plot_abund <- function(abund_summary, pdf_file, png_file) {

  #abund_summary <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/04_annotation/diamond/diamond_blastx_abundance.bac"

  summary_df = read.csv(abund_summary, sep="\t", header=T)

  species_summary <- summary_df %>%
    group_by(Sample, Species) %>%
    summarise(Sum_Mapped = sum(Reads_Mapped)) %>%
    top_n(10, Sum_Mapped)
  
  species_summary <- subset(species_summary, Species != "<not present>")
  
  # if a species is not in a particular sample, ggplot draws a wide column
  # want to have a 0 in this case so the column width is constant
  species_summary_wide = spread(species_summary, Sample, Sum_Mapped, fill=0)
  species_summary_long =  gather(species_summary_wide, Sample, Sum_Mapped, 2:ncol(species_summary_wide))

    # plot the summary table
  p <- ggplot(species_summary_long, aes(x=reorder(Species, Sum_Mapped), y=Sum_Mapped, fill=Sample)) +
    geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    scale_y_continuous(labels = comma) +
    ylab("Reads mapped") +
    coord_flip()

  # dynamically change figure height depending on number of samples
  # add 2 inches for every additional sample
  # and 0.1 inch for every speceis (sometimes less than 10 are recorded)
  spp = 0.1 * length(unique(species_summary$Species))
  ht = spp + (2 * (length(unique(summary_df$Sample))))
  
  ggsave(pdf_file, p, width=8, height=ht)
  ggsave(png_file, width=8, height=ht, dpi=300)

}

args <- commandArgs(trailingOnly = TRUE)
abund_summary = args[1]
pdf_file = args[2]
png_file = args[3]

plot_abund(abund_summary, pdf_file, png_file)


