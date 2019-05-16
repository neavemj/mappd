# plot the time component in the summary benchmarks file
# produced by 'summarise_benchmarks.py'

library(ggplot2)


plot_bench_time <- function(bench_summary, pdf_file, png_file) {
  
  #bench_summary <- "/datasets/work/AAHL_PDNGS_WORK/test_data/freshwater_prawn/benchmarks/tmp_benchmarks/tmp.txt"
  
  bench_fl = read.table(bench_summary, sep="\t", header=T)
  
  # not summarising like this at the moment 
  ## this bit just gets total time values for all samples, and all processes ##
  #total_time_process <- aggregate(bench_fl$s, by=list(sum=bench_fl$module), FUN=sum)
  #total_time_process$sample <- "combined"

  #total_time_sample <- aggregate(bench_fl$s, by=list(sample=bench_fl$sample), FUN=sum)
  #total_time_sample$sum <- "total_sample_time"

  #grand_total <- aggregate(total_time_sample$x, by=list(sum=total_time_sample$sum), FUN=sum)
  #grand_total$sample <- "combined"

  #total_df <- rbind(total_time_process, total_time_sample, grand_total) 
  ## ##
  
  ggplot(bench_fl, aes(x=s, y=process, color=sample)) +
    geom_point(size=2) +
    #geom_point(data=total_df, aes(x=x, y=sum), size=2) +
    scale_x_time() +
    xlab("time taken") +
    facet_grid(module ~ ., scales="free", space="free") +
    theme(strip.text.y = element_text(angle=0))

  # plot 1 inch per process
  fig_height = 1 * (length(unique(bench_fl$process)))
  
  ggsave(pdf_file, width=8, height=fig_height)
  ggsave(png_file, width=8, height=fig_height, dpi=300)
  
}


args <- commandArgs(trailingOnly = TRUE)
bench_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_bench_time(bench_summary, pdf_file, png_file)


