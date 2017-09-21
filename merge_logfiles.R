#!/usr/bin/env Rscript

# Read in package to read in command-line options.
library("optparse")

version <- "1.0"

option_list <- list(
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Comma-delimited list of logfiles to combine (required)." , 
              metavar="path"),
  
  make_option(c("-d", "--delim"), type="character", default="\t",
              help="Character specifying how logfiles are delimited (default: \"\t\").", 
              metavar="path"),
  
  make_option(c("-n", "--names"), type="character", default=NULL,
              help=paste("Optional comma-delimited strings that should be added to distinguish",
                         "the columns by file (default=NULL)."), 
              metavar="path"),
  
  make_option(c("-o", "--output"), type="character", default="combined_log.txt",
              help="Path to output file (default: \"combined_log.txt\").", 
              metavar="path"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
  
)

opt_parser <- OptionParser(option_list=option_list, 
                           usage = "%prog [options] -i log1.txt,log2.txt",
                           
                           description = paste("Basic script to read in multiple logfiles and to merge them by rows.\n",
                                               "The first column of each file is assumed to contain the sample names", sep="")
)

opt <- parse_args(opt_parser)

# Print out version if --version flag set.
if (opt$version) {
  cat("Wrapper version:", version, "\n")
  options_tmp <- options(show.error.messages=FALSE) 
  on.exit(options(options_tmp)) 
  stop()
}

if(is.null(opt$input)) {
  stop("paths to input logfiles need to be set.")
}

in_files <- strsplit(opt$input, ",")

# Check if only one file given.
if(length(in_files[[1]]) == 1) {
  stop("only one input file given.") 
}

# Check if names option given and make sure same length as files if so.
if(! is.null(opt$names)){
  in_names <- strsplit(opt$names, ",")
  
  if(length(in_files[[1]]) != length(in_names[[1]])) {
    stop("different numbers of infiles and names given.")
  }
}

# Define function to add new empty rows and give them specified rownames in df.
add_NA_rows <- function(rownames2add, df_in) {
  
  rows2add <- data.frame(matrix(NA, ncol=ncol(df_in), nrow=length(rownames2add)))
  colnames(rows2add) <- colnames(df_in)
  rownames(rows2add) <- rownames2add
  
  return(rbind(rows2add, df_in))
  
}

# Define function to merge 2 dataframes on rownames. 
# Will set any rows found in only 1 df to be all NAs in the df where it is missing.
cbind_w_missing <- function(df1, df2) {
 
   df1_only <- rownames(df1)[which(! rownames(df1) %in% rownames(df2))]
   df2_only <- rownames(df2)[which(! rownames(df2) %in% rownames(df1))]
  
   if(length(df1_only >= 1)){
     df2 <- add_NA_rows(df1_only, df2)
   }
   
   if(length(df2_only >= 1)){
     df1 <- add_NA_rows(df2_only, df1)
   }
   
   return(cbind(df1, df2))
   
}

# Keep file count marker.
file_count = 0

# Loop over index of all infiles.
for(i in 1:length(in_files[[1]])){
  infile <- read.table(in_files[[1]][i], header=TRUE, quote="", row.names=1, 
                       sep=opt$delim, stringsAsFactors=FALSE)
  
  if (! is.null(opt$names)) {
    colnames(infile) <- paste(in_names[[1]][i], colnames(infile), sep=".")
  }
  
  file_count = file_count + 1
  
  if(file_count == 1) {
    combined <- infile
  } else {
    combined <- cbind_w_missing(combined, infile)
  }
}

all_col <- colnames(combined)
combined$sample <- rownames(combined)
combined <- combined[,c("sample", all_col)]

write.table(x = combined, file = opt$output, quote = FALSE, sep = opt$delim,
            col.names = TRUE, row.names = FALSE)
