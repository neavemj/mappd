
# script to plot the overall results for trimming, rRNA mapping
# and euk, bac and vir abundances

library(ggplot2)
library(tidyverse)
library(scales)


plot_overall <- function(trim, rRNA, abund, pdf_file, png_file) {


  #trim <- "/datasets/work/AAHL_PDNGS_WORK/test_data/evag_RRV/logs/trimmomatic_PE/trim_logs.summary"
  #rRNA <- "/datasets/work/AAHL_PDNGS_WORK/test_data/evag_RRV/logs/rRNA_mapping_summary.tsv"
  #abund <- "/datasets/work/AAHL_PDNGS_WORK/test_data/evag_RRV/04_annotation/diamond/diamond_blastx_abundance.all"

  trim_df = read.csv(trim, sep="\t", header=T)
  # I'm not using reads if their mate is discarded
  # will combine these results into the 'dropped' category and convert to reads, not pairs
  trim_df$low_quality <- (trim_df$input_pairs - trim_df$both_surviving) * 2
  colnames(trim_df) <- c("Sample", "input_pairs", "both_surviving", "forward_only", "reverse_only", "dropped", "low_quality")

  # make long format
  trim_long =  gather(trim_df, Type, Reads, low_quality)
  trim_long <- trim_long[,c("Sample", "Type", "Reads")]

  rRNA_df = read.csv(rRNA, sep="\t", header=F, row.names=NULL)
  # make headers match for later rbind
  colnames(rRNA_df) <- c("Sample", "Type", "Reads")
  rRNA_df <- subset(rRNA_df, Type!="mRNA_reads")

  abund_df = read.csv(abund, sep="\t", header=T, row.names=NULL)

  euk_df = subset(abund_df, Kingdom=="Eukaryota")
  bac_df = subset(abund_df, Kingdom=="Bacteria")
  vir_df = subset(abund_df, Kingdom=="Viruses")

  # function to check that taxa dataframes have something in them
  # return a dataframe with 0 reads if not
  process_taxa <- function(taxa_df, name){
    if(length(taxa_df$Sample) > 0){
      summary <- taxa_df %>%
        group_by(Sample) %>%
        summarise(Reads = sum(Reads_Mapped))
      summary["Type"] <- name
    } else {
      num = length(unique(abund_df$Sample))
      summary <- data.frame(Sample <- unique(abund_df$Sample),
                            Reads <- rep(0, num),
                            Type <- rep(name, num))
      colnames(summary) <- c("Sample", "Reads", "Type")
    }
    return(summary)
  }

  euk_summary <- process_taxa(euk_df, "Eukaryote")
  bac_summary <- process_taxa(bac_df, "Bacteria")
  vir_summary <- process_taxa(vir_df, "Virus")

  overall_df <- rbind(trim_long, rRNA_df, euk_summary, bac_summary, vir_summary)

  # figure out the total number of annotated / unannotated reads from these numbers
  # should verify this by looking at:
    # 1) reads that didn't map to the contigs (ie. didn't form contigs)
    # 2) plus reads that mapped to contigs that didn't get annotated
    # the sum of these two things should equal the Unannotated calcs below
  unannot_df <- overall_df %>%
    group_by(Sample) %>%
    summarise(Total_Annotated = sum(Reads))

  unannot_df <- merge(trim_df, unannot_df, by="Sample")

  unannot_df$Reads <- (unannot_df$input_pairs * 2) - unannot_df$Total_Annotated
  print(unannot_df)
  # if most of the reads have been annotated, this can sometimes be less than 0
  # due to paired end calculations in the bowtie output files
  # could improve this by using a SQL database and tracking individual reads
  #if(unannot_df$Reads < 0){
  #  unannot_df$Reads <- 0
  #}

  unannot_df$Type <- "Unannotated"
  unannot_df <- unannot_df[,c("Sample", "Type", "Reads")]

  overall_df <- rbind(overall_df, unannot_df)

  # make the colours a bit more sensible
  cols <- c("low_quality" = "#a6761d", "rRNA_LSU" = "#FDBF6F", "rRNA_SSU" = "#FF7F00", "Eukaryote" = "#1b9e77", "Bacteria" = "#7570b3", "Virus" = "#e7298a", "Unannotated" = "#666666")

  # order the categories
  # flip everything around due to my coord_flip() in ggplot call

  overall_df$Type <- factor(overall_df$Type, levels=rev(c("Virus", "Bacteria", "Eukaryote", "rRNA_LSU", "rRNA_SSU", "low_quality", "Unannotated")))

  overall_df$Sample <- as.factor(overall_df$Sample)
  overall_df$Sample <- factor(overall_df$Sample, levels=rev(levels(overall_df$Sample)))

    # plot the summary table
  p <- ggplot(overall_df, aes(x=Sample, y=Reads, fill=Type)) +
    geom_bar(stat='identity') +
    theme_bw() +
    scale_fill_manual(values=cols) +
    theme(axis.title.y = element_blank()) +
    scale_y_continuous(labels = comma) +
    ylab("Reads") +
    coord_flip() +
    guides(fill = guide_legend(reverse=T))

  # dynamically change figure height depending on number of samples
  # although, it has to be at least 3 inches high for the legend
  # add 1 inch for every additional sample

  if(length(unique(overall_df$Sample)) < 3){
    ht = 3
  } else {
    ht = 1 * (length(unique(overall_df$Sample)))
  }

  ggsave(pdf_file, p, width=8, height=ht)
  ggsave(png_file, width=8, height=ht, dpi=300)
  write.table(overall_df, "overall_counts.txt", sep="\t", quote=F, row.names=F)

}

args <- commandArgs(trailingOnly = TRUE)
trim = args[1]
rRNA = args[2]
abund = args[3]
pdf_file = args[4]
png_file = args[5]

plot_overall(trim, rRNA, abund, pdf_file, png_file)



