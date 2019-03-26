# plot the memory component from the summary benchmarks file
# produced by 'summarise_benchmarks.py'

library(ggplot2)


plot_bench_mem <- function(bench_summary, pdf_file, png_file) {
  
  bench_fl = read.table(bench_summary, sep="\t", header=T)

  ggplot(bench_fl, aes(x=max_vms, y=process, color=sample)) +
    geom_point(size=2) +
    expand_limits(x=0)

  fig_height = 1.5 + (length(unique(bench_fl$process)) / 3)
  
  ggsave(pdf_file, width=8, height=fig_height)
  ggsave(png_file, width=8, height=fig_height, dpi=300)
  
}


args <- commandArgs(trailingOnly = TRUE)
bench_summary = args[1]
pdf_file = args[2]
png_file = args[3]
plot_bench_mem(bench_summary, pdf_file, png_file)


