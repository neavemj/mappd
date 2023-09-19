
library(rmarkdown)
library(optparse)

# use optparse to grab command line arguments for the report

option_list <- list(
  # required args
  make_option(c("--config_file"), type="character", default=NULL,
              help="Config file used to run the pipeline", metavar="character"),
              
  make_option(c("--dag"), type="character", default=NULL,
              help="PNG DAG graph of the ran pipeline", metavar="character"),

  make_option(c("--software_list"), type="character", default=NULL,
              help="all software and versions from conda list", metavar="character"),

  make_option(c("--bench_time"), type="character", default=NULL,
              help="PNG graph of time benchmarks", metavar="character"),

  make_option(c("--bench_mem"), type="character", default=NULL,
              help="PNG graph of memory benchmarks", metavar="character"),

  make_option(c("--tech_summary"), type="character", default=NULL,
              help="table of technical summary", metavar="character"),

  make_option(c("--trim"), type="character", default=NULL,
              help="logs of the trim summary", metavar="character"),

  make_option(c("--rRNA"), type="character", default=NULL,
              help="logs of the rRNA mapping summary", metavar="character"),

  make_option(c("--abund"), type="character", default=NULL,
              help="the diamond_blastx_abundance.all file", metavar="character"),

  make_option(c("--taxa_figures"), type="character", default=NULL,
              help="file containing which taxa figures were produced", metavar="character"),

  make_option(c("--sample_abundances"), type="character", default=NULL,
              help="taxa abundance tables for each sample", metavar="character"),

  make_option(c("--output"), type="character", default=NULL,
              help="name and file path for html report", metavar="character"),

  make_option(c("--output_dir"), type="character", default=NULL,
              help="directory to put html file. Otherwise rmarkdown puts
              the output in the same directory as the markdown document", metavar="character"),
             
  make_option(c("--rmarkdown"), type="character", default=NULL,
              help="location of rmarkdown file", metavar="character")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# need to grab the working project directory here
# otherwise rmarkdown uses the directory where the markdown document
# is located as the working directory. 

working_dir <- getwd()

# now read data for making the report

render(opt$rmarkdown,
    params = list(
    config_file = opt$config_file,
    dag = opt$dag,
    software_list = opt$software_list,
    bench_time = opt$bench_time,
    bench_mem = opt$bench_mem,
    tech_summary = opt$tech_summary,
    trim = opt$trim,
    rRNA = opt$rRNA,
    abund = opt$abund,
    taxa_figures = opt$taxa_figures,
    sample_abundances = opt$sample_abundances
    ),
    output_file = opt$output,
    output_dir = opt$output_dir)

