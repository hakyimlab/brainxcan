library(optparse)

option_list <- list(
    make_option(c("-i", "--input_prefix"), type="character", default=NULL,
                help="The prefix of all these BrainXcan analysis results",
                metavar="character"),
    # make_option(c("-m", "--idp_meta_file"), type="character", default=NULL,
    #             help="A meta file for annotating IDPs",
    #             metavar="character"),
    make_option(c("-c", "--color_code_yaml"), type="character", default=NULL,
                help="Color coding",
                metavar="character"),
    make_option(c("-r", "--rlib"), type="character", default=NULL,
                help="The path to report helper functions",
                metavar="character"),
    make_option(c("-n", "--ntop"), type="numeric", default=NULL,
                help="Number of top IDP associations to show",
                metavar="character"),
    make_option(c("-p", "--phenotype_name"), type="character", default=NULL,
                help="Phenotype name to show in the report",
                metavar="character"),
    make_option(c("-t", "--rmd_template"), type="character", default=NULL,
                help="R Markdown template",
                metavar="character"),
    make_option(c("-o", "--output_html"), type="character", default=NULL,
                help="Output HTML filename",
                metavar="character"),
    make_option(c("-k", "--bxcan_pval_col"), type="character", default=NULL,
                help="The name of the p-value column to be used in BrainXcan results.",
                metavar="character"),
    make_option(c("--fdr_cutoff"), type="numeric", default=NULL,
                help="FDR cutoff to call significant",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

params = list(
  color_code_yaml = opt$color_code_yaml,
  input_prefix = opt$input_prefix,
  rlib = opt$rlib,
  phenotype_name = opt$phenotype_name,
  ntop = opt$ntop,
  bxcan_pval_col = opt$bxcan_pval_col,
  fdr_cutoff = opt$fdr_cutoff
)

rmd_file <- paste0(tools::file_path_sans_ext(opt$output_html), '.Rmd')
system(paste('cp', opt$rmd_template, rmd_file))

rmarkdown::render(
  rmd_file, 
  params = params, 
  envir = new.env(),
  output_dir = dirname(opt$output_html),
  output_file = basename(opt$output_html)
)
