
# script to plot the overall results for trimming, rRNA / host mapping
# and euk, bac and vir abundances

library(ggplot2)
library(tidyverse)
library(scales)


plot_overall <- function(trim, rRNA_host, euk, bac, vir, pdf_file, png_file) {

  
  trim <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/logs/trimmomatic_PE/trim_logs.summary"
  rRNA_host <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/logs/mapping_summary.tsv"
  euk <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/04_annotation/diamond/diamond_blastx_abundance.euk"
  bac <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/04_annotation/diamond/diamond_blastx_abundance.bac"
  vir <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/04_annotation/diamond/diamond_blastx_abundance.vir"

  trim_df = read.csv(trim, sep="\t", header=T)
  # I'm not using reads if they're mate is discarded
  # will combine these results into the 'dropped' category and convert to reads, not pairs
  trim_df$low_quality <- (trim_df$input_pairs - trim_df$both_surviving) * 2
  colnames(trim_df) <- c("Sample", "input_pairs", "both_surviving", "forward_only", "reverse_only", "dropped", "low_quality")
  
  # make long format
  trim_long =  gather(trim_df, Type, Reads, low_quality)
  trim_long <- trim_long[,c("sample", "Type", "Reads")]
  colnames(trim_long) <- c("Sample", "Type", "Reads")
  
  rRNA_df = read.csv(rRNA_host, sep="\t", header=T)
  # make headers match for later rbind
  rRNA_df$Reads <- rRNA_df$Paired_Reads * 2
  rRNA_df <- rRNA_df[,c("Sample", "Type", "Reads")]
  # the 'host' here should actually be combined with the Eukaryotes
  # should make sure that my host module always uses a Euk
  # or aborts otherwise
  levels(rRNA_df$Type)[levels(rRNA_df$Type)=="host"] <- "Eukaryotic"
  # don't need the 'surviving' number
  rRNA_df <- subset(rRNA_df, Type!="surviving")
  
  
  euk_df = read.csv(euk, sep="\t", header=T)
  bac_df = read.csv(bac, sep="\t", header=T)
  vir_df = read.csv(vir, sep="\t", header=T)

  euk_summary <- euk_df %>%
    group_by(Sample) %>%
    summarise(Reads = sum(Reads_Mapped))
  euk_summary["Type"] <- "Eukaryotic"

  
  bac_summary <- bac_df %>%
    group_by(Sample) %>%
    summarise(Reads = sum(Reads_Mapped))
  bac_summary["Type"] <- "Bacteria"
  
  vir_summary <- vir_df %>%
    group_by(Sample) %>%
    summarise(Reads = sum(Reads_Mapped))
  vir_summary["Type"] <- "Viruses"
  
  
  overall_df <- rbind(trim_long, rRNA_df, euk_summary, bac_summary, vir_summary)

  unannot_df <- overall_df %>%
    group_by(Sample) %>%
    summarise(Total_Annotated = sum(Reads))

  unannot_df <- merge(trim_df, unannot_df, by="Sample")
  unannot_df$Reads <- (unannot_df$input_pairs * 2) - unannot_df$Total_Annotated
  unannot_df$Type <- "Unannotated"
  unannot_df <- unannot_df[,c("Sample", "Type", "Reads")]
  
  overall_df <- rbind(overall_df, unannot_df)
  
    # plot the summary table
  p <- ggplot(overall_df, aes(x=Sample, y=Reads, fill=Type)) +
    geom_bar(stat='identity') +
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


