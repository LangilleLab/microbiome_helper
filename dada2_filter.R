#!/usr/bin/env Rscript

# Read in package to read in command-line options.
library("optparse")

version <- "1.0"

option_list <- list(

  make_option(c("-f", "--f_path"), type="character", default=NULL,
              help="Location of forward reads (required)." , metavar="path"),
  
  make_option(c("-r", "--r_path"), type="character", default=NULL,
              help="Location of reverse reads (if applicable). Same as forward reads by default.",
              metavar="path"),
  
  make_option(c("-s", "--single"), action = "store_true", type="logical", default=FALSE,
              help="Flag to indicate that reads are single-end (default: FALSE).", metavar = "boolean"),
  
  make_option(c("--f_match"), type="character", default="_R1_.*fastq.*",
              help="Regular expression to match in forward reads' filenames (default: _R1_.*fastq.*).", 
              metavar="PATTERN"),
  
  make_option(c("--r_match"), type="character", default="_R2_.*fastq.*",
              help="Regular expression to match in reverse reads' filenames (default: \"_R2_.*fastq.*\").",
              metavar="PATTERN"),
  
  make_option(c("--sample_delim"), type="character", default="_",
              help=paste("Character to split filenames on to determine sample name (default: \"_\").",
              "Sample names are assumed to be field 1 after splitting.", sep=" "),
              metavar="PATTERN"),
  
  make_option(c("--truncQ"), type="character", default="2",
              help="Reads will be truncated at the first instance of Q <= truncQ (default: 2).", 
              metavar="numeric,numeric"),
  
  make_option(c("--truncLen"), type="character", default="0",
              help="Length at which to truncate reads - shorter reads will be discarded (default: 0).", 
              metavar="numeric,numeric"),
  
  make_option(c("--trimLeft"), type="character", default="0",
              help="Number of nucleotides to trim from start of reads (default: 0).", 
              metavar="numeric,numeric"),
  
  make_option(c("--maxLen"), type="character", default="Inf",
              help="Maximum read length BEFORE trimming and truncation (default: Inf).", 
              metavar="numeric,numeric"),
  
  make_option(c("--minLen"), type="character", default="20",
              help="Minimum read length AFTER trimming and truncation (default: 20).", 
              metavar="numeric,numeric"),
  
  make_option(c("--maxN"), type="character", default="0",
              help="Maximum number of Ns after truncation - should be 0 for use with DADA2 (default: 0).", 
              metavar="numeric,numeric"),

  make_option(c("--minQ"), type="character", default="0",
              help="Minimum quality score allowed in read after truncation (default: 0).", 
              metavar="numeric,numeric"),
  
  make_option(c("--maxEE"), type="character", default="Inf",
              help="Maximum expected errors allowed in a read after truncation (default: Inf). EE = sum(10^(-Q/10)).", 
              metavar="numeric,numeric"),  
  
  make_option(c("-t", "--threads"), type="integer", default=1,
              help="Number of threads to use (default: 1).", metavar="integer"),
  
  make_option(c("-n", "--num_reads"), type="integer", default=1e5,
              help="Number of reads to read into memory at once (default: 1e5).", metavar="integer"),  
  
  make_option(c("-o", "--output"), type="character", default="filtered_fastqs",
              help="Output folder for filtered reads (default: filtered_fastqs).", metavar="path"),
  
  make_option(c("--log"), type="character", default="dada2_filter_read_counts.txt",
              help="Output logfile for read count table (default: dada2_filter_read_counts.txt).", metavar="path"),
  
  make_option(c("--no_rm_phiX"), action = "store_true", type="logical", default=FALSE,
              help="Flag to indicate that reads matching PhiX genome should NOT be discarded (default: FALSE).", 
              metavar = "boolean"),
  
  make_option(c("--no_gzip"), action = "store_true", type="logical", default=FALSE,
              help="Flag to indicate that output FASTQs should NOT be gzipped (default: FALSE).", 
              metavar = "boolean"),
  
  make_option(c("--matchIDs"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that paired-end reads that don't share ids",
                         "in both the forward and reverse FASTQs will be excluded (default: FALSE)", sep=" "), 
              metavar = "boolean"),
  
  make_option(c("--id_sep"), type="character", default="\\s",
              help="The separator between fields in the id line for paired-end FASTQs (default: \"\\s\").", 
              metavar = "character"),
  
  make_option(c("--id_field"), type="character", default=NULL,
              help=paste("Field of the FASTQ id line containing sequenced id. By default",
                         "this will be automatically detected, which will work for most users.",
                         "Only for paired-end data.", sep=" "), 
              metavar = "boolean"),
  
  make_option(c("--verbose"), action = "store_true", type="logical", default=FALSE,
              help="Write out status messages (default: FALSE).", metavar = "boolean"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")

)

opt_parser <- OptionParser(option_list=option_list, 
                           usage = "%prog [options] -f PATH",
                           
                           description = paste("\nThis is a wrapper for DADA2 that is based on the authors\'", 
                                               "big data tutorial available here:",
                                               "https://benjjneb.github.io/dada2/bigdata.html.\n\nBe sure to cite the DADA2",
                                               "paper if you use this script:\nCallahan BJ et al. 2016. DADA2:",
                                               "High-resolution sample inference from Illumina amplicon data.",
                                               "Nature Methods 13:581-583.\n\nNote this script was tested with",
                                               "DADA2 v1.4.0\n\nUSAGE NOTE: when passing filtering parameters for paired-end", 
                                               "reads if you pass one value it will be used for both forward and reverse reads.",
                                               "\nIf you pass two values separated by a comma then these values will be for the forward",
                                               "and reverse reads respectively.", sep=" ")
                          )

opt <- parse_args(opt_parser)

# Function to parse DADA2 filtering params and to throw error if any unexpected input.
parse_dada2_filt_params <- function(in_param, single_end, param_name) {
  
  # Split by comma (in case separate options given for R1 and R2).
  in_param_split <- as.vector(strsplit(in_param, ",")[[1]])
  
  # Report error if empty string.
  if(length(in_param_split) == 0) {
    stop(paste("Error when parsing parameter", param_name, "-",
               "no values detected.", sep=" "))
  }
  
  # Make sure that only one value given if SE.
  if(single_end & length(in_param_split) != 1) {
        stop(paste("Error when parsing parameter", param_name, "-",
                   "expected 1 input value, but got", length(in_param_split), 
                   sep=" "))
             
  # Make sure that max of 2 values given if PE.
  } else if(! single_end & length(in_param_split) > 2) {
    stop(paste("Error when parsing parameter", param_name, "-",
               "expected max of 2 input values, but got", length(in_param_split), 
               sep=" "))
  }
  
  if (length(in_param_split) == 1) {
    return(as.numeric(in_param_split[1]))
  } else {
    return(as.numeric(in_param_split[1]), as.numeric(in_param_split[2]))
  }
  
}



# Print out version if --version flag set.
if (opt$version) {
  cat("Wrapper version:", version, "\n")
  options_tmp <- options(show.error.messages=FALSE) 
  on.exit(options(options_tmp)) 
  stop()
}

# Check if path to forward files set, if not stop job.
if(is.null(opt$f_path)) {
  stop("path to forward FASTQs needs to be set.")
}

# Get forward and reverse filenames.
forward_in <- sort(list.files(opt$f_path, pattern=opt$f_match, full.names=TRUE))
forward_samples <- sapply(strsplit(basename(forward_in), opt$sample_delim), `[`, 1)
forward_out <- file.path(opt$output, paste0(forward_samples, "_R1_filt.fastq.gz"))

# Convert input filtering params to vectors and check for obvious input errors.
filt_params_raw <- opt[c("truncQ", "truncLen", "trimLeft", "maxLen",
                         "minLen", "maxN", "minQ", "maxEE")]
filt_params <- list(c())
for (param in names(filt_params_raw)) {
  filt_params[[param]] <- parse_dada2_filt_params(in_param=filt_params_raw[[param]],
                                                  single_end=opt$single,
                                                  param_name=param)
}

# Set multithread option (FALSE if no, otherwise give # core).
if (opt$threads > 1) {
  multithread_opt <- opt$threads
} else {
  multithread_opt <- FALSE
}

# Finally read in and print out DADA2 version.
library("dada2")

if(opt$verbose) {
  cat("Running DADA2 version:", as.character(packageVersion("dada2")), "\n")
  cat("Running filterAndTrim function with the below options.\n")
}

# Set id_field output log to be character of "NULL" so that it is printed (if applicable).
id_field_log <- opt$id_field
if(is.null(opt$id_field)) {
  id_field_log <- "NULL"
}

# Run DADA2 on just forward reads in single-end data.
if (opt$single) {
   
  if(opt$verbose) {
    
    cat("Single-end mode\n\n")
    cat(paste("fwd=", paste(forward_in, collapse=",") , "\n", 
              "filt=", paste(forward_out, collapse=","), "\n",
              "compress=", !opt$no_gzip, "\n",
              "truncQ=", filt_params$truncQ, "\n",  
              "truncLen=", filt_params$truncLen, "\n", 
              "trimLef=", filt_params$trimLeft, "\n",  
              "maxLen=", filt_params$maxLen, "\n", 
              "minLen=", filt_params$minLen, "\n",  
              "maxN=", filt_params$maxN, "\n", 
              "minQ=", filt_params$minQ, "\n",
              "maxEE=", filt_params$maxEE, "\n",
              "rm.phix=", !opt$no_rm_phiX, "\n",
              "multithread=", multithread_opt, "\n",
              "n=", opt$num_reads, "\n",
              "verbose=", opt$verbose, "\n\n", sep=""))
  }
  
  read_counts <- filterAndTrim(fwd=forward_in, filt=forward_out, compress=!opt$no_gzip,
                               truncQ=filt_params$truncQ, truncLen=filt_params$truncLen,
                               trimLeft=filt_params$trimLeft, maxLen=filt_params$maxLen,
                               minLen=filt_params$minLen, maxN=filt_params$maxN,
                               minQ=filt_params$minQ, maxEE=filt_params$maxEE, 
                               rm.phix=!opt$no_rm_phiX, multithread=multithread_opt,
                               n=opt$num_reads, verbose=opt$verbose)

} else {
  
  # Read in path to reverse FASTQs.
  if(is.null(opt$r_path)) {
    opt$r_path <- opt$f_path
  }
  
  reverse_in <- sort(list.files(opt$r_path, pattern=opt$r_match, full.names=TRUE))
  reverse_samples <- sapply(strsplit(basename(reverse_in), opt$sample_delim), `[`, 1)
  reverse_out <- file.path(opt$output, paste0(reverse_samples, "_R2_filt.fastq.gz"))
  
  # Check that sample names for forward and reverse FASTQs are the same.
  if(! identical(forward_samples, reverse_samples)){
    stop(paste("Sample names parsed from forward and reverse filenames don't match.",
               "\n\nForward sample names:", forward_samples,
               "\n\nReverse sample names:", reverse_samples))
  }
  
  if(opt$verbose) {
    cat("Paired-end mode\n")
    cat(paste("fwd=", paste(forward_in, collapse=",") , "\n", 
              "filt=", paste(forward_out, collapse=","), "\n",
              "rev=", paste(reverse_in, collapse=","), "\n",
              "filt.rev=", paste(reverse_out, collapse=","), "\n",
              "compress=", !opt$no_gzip, "\n",
              "truncQ=", filt_params$truncQ, "\n",  
              "truncLen=", filt_params$truncLen, "\n", 
              "trimLef=", filt_params$trimLeft, "\n",  
              "maxLen=", filt_params$maxLen, "\n", 
              "minLen=", filt_params$minLen, "\n",  
              "maxN=", filt_params$maxN, "\n", 
              "minQ=", filt_params$minQ, "\n",
              "maxEE=", filt_params$maxEE, "\n",
              "rm.phix=", !opt$no_rm_phiX, "\n",
              "multithread=", multithread_opt, "\n",
              "n=", opt$num_reads, "\n",
              "matchIDs=", opt$matchIDs, "\n", 
              "id.sep=", opt$id_sep, "\n",
              "id.field=", opt$id_field, "\n", 
              "verbose=", opt$verbose, "\n\n", sep=""))
  }
  
  read_counts <- filterAndTrim(fwd=forward_in, filt=forward_out, rev=reverse_in,
                               filt.rev=reverse_out, compress=!opt$no_gzip,
                               truncQ=filt_params$truncQ, truncLen=filt_params$truncLen,
                               trimLeft=filt_params$trimLeft, maxLen=filt_params$maxLen,
                               minLen=filt_params$minLen, maxN=filt_params$maxN,
                               minQ=filt_params$minQ, maxEE=filt_params$maxEE, 
                               rm.phix=!opt$no_rm_phiX, multithread=multithread_opt,
                               n=opt$num_reads, matchIDs=opt$matchIDs, id.sep=opt$id_sep,
                               id.field=opt$id_field, verbose=opt$verbose)
  
}

# Convert log table rownames from R1 FASTQs to sample names.
read_counts <- as.data.frame(read_counts)
rownames(read_counts) <- sapply(strsplit(rownames(read_counts), opt$sample_delim), `[`, 1)
colnames(read_counts) <- c("input", "filtered")

# Print out input and output read counts to logfile.
write.table(x = read_counts, file = opt$log, quote = FALSE, sep="\t",
            col.names = NA, row.names = TRUE)
