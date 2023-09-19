
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

  make_option(c("--output"), type="character", default=NULL,
              help="name and file path for html report", metavar="character"),

  make_option(c("--output_dir"), type="character", default=NULL,
              help="directory to put html file. Otherwise rmarkdown puts
              the output in the same directory as the markdown document", metavar="character"),
             
  make_option(c("--rmarkdown"), type="character", default=NULL,
              help="location of rmarkdown file", metavar="character")
)


opt_parser <- OptionParser(option_list=option_list)
# need to specify positional_arguments as true because of the sample abundances files
# this allows me to pass a variable number of multiple files (don't know how many samples beforehand)
opt <- parse_args(opt_parser, positional_arguments=TRUE)

# the normal options go into a opt$options variable
# positional arguments (i.e. sample abundance files) go into opt$args
#print(opt$options)
#print(opt$args)

# need to grab the working project directory here
# otherwise rmarkdown uses the directory where the markdown document
# is located as the working directory. 

working_dir <- getwd()

# now read data for making the report

render(opt$options$rmarkdown,
    params = list(
    config_file = opt$options$config_file,
    dag = opt$options$dag,
    software_list = opt$options$software_list,
    bench_time = opt$options$bench_time,
    bench_mem = opt$options$bench_mem,
    tech_summary = opt$options$tech_summary,
    trim = opt$options$trim,
    rRNA = opt$options$rRNA,
    abund = opt$options$abund,
    taxa_figures = opt$options$taxa_figures,
    sample_abundances = opt$args
    ),
    output_file = opt$options$output,
    output_dir = opt$options$output_dir)

