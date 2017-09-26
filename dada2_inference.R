#!/usr/bin/env Rscript

# Read in package to read in command-line options.
library("optparse")

version <- "1.0"

option_list <- list(

  make_option(c("-f", "--f_path"), type="character", default=NULL,
              help="Location of forward reads (required)." , metavar="path"),
  
  make_option(c("--seed"), type="integer", default=NULL,
              help="Random seed to make command reproducible.", metavar = "integer"),
  
  make_option(c("-r", "--r_path"), type="character", default=NULL,
              help="Location of reverse reads (if applicable). Same as forward reads by default.",
              metavar="path"),
  
  make_option(c("-s", "--single"), action = "store_true", type="logical", default=FALSE,
              help="Flag to indicate that reads are single-end (default: FALSE).", metavar = "boolean"),

  make_option(c("--f_match"), type="character", default="_R1_.*fastq.*",
              help="Regular expression to match in forward reads' filenames (default: \"_R1_.*fastq.*\").", 
              metavar="PATTERN"),
  
  make_option(c("--r_match"), type="character", default="_R2_.*fastq.*",
              help="Regular expression to match in reverse reads' filenames (default: \"_R2_.*fastq.*\").",
              metavar="PATTERN"),
  
  make_option(c("--sample_delim"), type="character", default="_",
              help=paste("Character to split filenames on to determine sample name (default: \"_\").",
              "Sample names are assumed to be field 1 after splitting.", sep=" "),
              metavar="PATTERN"),
  
  make_option(c("-t", "--threads"), type="integer", default=1,
              help="Number of threads to use (default: 1).", metavar="integer"),
  
  make_option(c("-o", "--output"), type="character", default="seqtab.rds",
              help="Output RDS file for sequence table (default: seqtab.rds).", metavar="path"),
  
  make_option(c("-n", "--num_learn"), type="integer", default=1e6,
              help=paste("Min. number of reads for error rate learning (default: 1e6).", 
                         "Samples will be read in until this read count is reached or all",
                         "samples are read in.", sep=" "), metavar="integer"),
  
  make_option(c("--num_derep"), type="integer", default=1e6,
              help="Number of reads to read into memory when dereplicating (default: 1e6).", 
                    metavar="integer"),
  
  make_option(c("--randomize"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that samples should be read in random order",
                         "for learnErrors step (default: FALSE).", sep=" "),
              metavar = "boolean"),
  
  make_option(c("--plot_errors"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that estimated error rates should be plotted",
                         "to pdfs (default: FALSE).", sep=" "),
              metavar = "boolean"),
  
  make_option(c("--selfConsist"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that dada algorithm should alternate between sample",
                         "inference and error rate estimate until convergence (default: FALSE).", 
                         sep=" "), metavar = "boolean"),
  
  make_option(c("--pool"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that all samples should be pooled before sample",
                         "inference step (default: FALSE).", sep=" "), metavar = "boolean"),
  
  make_option(c("--minOverlap"), type="integer", default=20,
              help="Min length of overlap required for merging reads (default: 20).", 
              metavar="integer"),
  
  make_option(c("--maxMismatch"), type="integer", default=0,
              help="Max number of mismatches allowed in overlap when merging reads (default: 0).", 
              metavar="integer"),
  
  make_option(c("--returnRejects"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that pairs that had more than the max number", 
                         "of mismatches should be retained in the output dataframe",
                         "(default: FALSE).", sep=" "), metavar = "boolean"),
  
  make_option(c("--justConcatenate"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that forward and reverse reads should just", 
                         "be concatenated with 10 Ns as spacers between them",
                         "(default: FALSE).", sep=" "), metavar = "boolean"),
  
  make_option(c("--trimOverhang"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that overhands in alignment between forward",
                         "and reverse reads are trimmed off. This can happen when the reads",
                         "are bigger than the amplicon (default: False).", sep=" "), 
              metavar = "boolean"),
  
  make_option(c("--log"), type="character", default="dada2_inferred_read_counts.txt",
              help="Output logfile for read count table (default: dada2_inferred_read_counts.txt).", 
              metavar="path"),
  
  make_option(c("--verbose"), action = "store_true", type="logical", default=FALSE,
              help="Write out status messages (default: FALSE).", metavar = "boolean"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")

)

opt_parser <- OptionParser(
                    option_list=option_list, 
                    usage = "%prog [options] -f PATH --seed 123",
                    description = paste(
                             "\nThis is a wrapper for the DADA2 inference step that is", 
                             "based on the authors\' big data tutorial available here:",
                             "https://benjjneb.github.io/dada2/bigdata.html.\n\nBe sure to cite the DADA2",
                             "paper if you use this script:\nCallahan BJ et al. 2016. DADA2:",
                             "High-resolution sample inference from Illumina amplicon data.",
                             "Nature Methods 13:581-583.\n\nNote this script was tested with",
                             "DADA2 v1.4.0", sep=" ")
                          )

opt <- parse_args(opt_parser)

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

# Set multithread option (FALSE if no, otherwise give # core).
if (opt$threads > 1) {
  multithread_opt <- opt$threads
} else {
  multithread_opt <- FALSE
}

# Get forward filenames.
forward_in <- sort(list.files(opt$f_path, pattern=opt$f_match, full.names=TRUE))
forward_samples <- sapply(strsplit(basename(forward_in), opt$sample_delim), `[`, 1)
names(forward_in) <- forward_samples

# Read in and print out DADA2 version.
library("dada2")

if(opt$verbose) {
  cat("Running DADA2 version:", as.character(packageVersion("dada2")), "\n")

}

# Set random seed for reproducibility if set.
if(! is.null(opt$seed)) {
  set.seed(opt$seed)
  if(opt$verbose){
    cat("Random seed:", as.character(opt$seed), "\n")
  }
}
  
if(opt$verbose) {
    cat("Running learnErrors with below options.\n")
    cat(paste("fls=", paste(forward_in, collapse=",") , "\n",
              "nreads=", opt$num_learn, "\n",
              "multithread=", multithread_opt, "\n",
              "randomize=", opt$randomize, "\n\n", sep=""))
}

# Learn error model for forward reads.
err_forward <- learnErrors(fls = forward_in, 
                   nreads = opt$num_learn, 
                   multithread = multithread_opt, 
                   randomize=opt$randomize)

# Initialize empty list with length of number of samples.
seq_variants <- vector("list", length(forward_samples))
names(seq_variants) <- forward_samples

if(opt$verbose) {
  cat("Will run derepFastq and dada with the below options.\n")
  cat("derepFastq n =", opt$num_derep, "\n")
  cat("dada selfConsist =", opt$selfConsist, "\n")
  cat("dada pool =", opt$pool, "\n")
  cat("dada multithread =", multithread_opt, "\n")
}

# Run dereplication and dada inference just on forward reads in single-end mode.
if(opt$single) {
  
  if(opt$verbose) {
    cat("Running dereplication and dada algorithm in single-end mode\n")
  }
  
  # Initialize dataframe to keep track of read counts.
  log_counts <- data.frame(matrix(NA, nrow=length(forward_samples), ncol=3))
  rownames(log_counts) <- forward_samples
  colnames(log_counts) <- c("sample", "derep_sum", "denoised")
  
  # Loop over all samples.
  for(sample in forward_samples) {

    if(opt$verbose) {
      cat("Processing: ", sample, "\n")
    }
    
    # Dereplicate reads.
    derep_forward <- derepFastq(forward_in[[sample]], 
                                n=opt$num_derep,
                                verbose=opt$verbose)
    
    # Run inference and add inferred sequences to output list.
    seq_variants[[sample]] <- dada(derep_forward, 
                                   err=err_forward, 
                                   multithread=multithread_opt)
    
    # Add read counts to log dataframe for debug purposes.
    log_counts[sample,] <- c(sample, sum(derep_forward$uniques), sum(seq_variants[[sample]]$denoised))

  }

# Learn errors for reverse too if paired.
} else {
  
  # Initialize dataframe to keep track of read counts.
  log_counts <- data.frame(matrix(NA, nrow=length(forward_samples), ncol=6))
  rownames(log_counts) <- forward_samples
  colnames(log_counts) <- c("sample", "forward_derep_sum", "reverse_derep_sum", 
                            "forward_denoised", "reverse_denoised", "merged")
  
  # Read in path to reverse FASTQs.
  if(is.null(opt$r_path)) {
    opt$r_path <- opt$f_path
  }
  
  # Read in reverse FASTQs
  reverse_in <- sort(list.files(opt$r_path, pattern=opt$r_match, full.names=TRUE))
  reverse_samples <- sapply(strsplit(basename(reverse_in), opt$sample_delim), `[`, 1)
  names(reverse_in) <- reverse_samples
  
  # Check if forward and reverse sample names match.
  if(! identical(forward_samples, reverse_samples)){
    stop(paste("\n\nSample names parsed from forward and reverse filenames don't match.",
               "\nForward sample name:", forward_samples,
               "\nReverse sample name:", reverse_samples,
               "\n\nUse the -s option if your reads are single-end.\n"))
  }
  
  if(opt$verbose) {
    cat("Running learnErrors with below options for reverse reads.\n",
        "fls =", paste(reverse_in, collapse=",") , "\n",
        "nreads =", opt$num_learn, "\n",
        "multithread =", multithread_opt, "\n",
        "randomize =", opt$randomize, "\n\n",
        "Running dereplication and dada algorithm in paired-end mode.\n\n",
        "Will run mergePairs with below options.\n",
        "minOverlap =", opt$minOverlap, "\n",
        "maxMismatch =", opt$maxMismatch, "\n",
        "returnRejects =", opt$returnRejects, "\n",
        "justConcatenate =", opt$justConcatenate, "\n",
        "trimOverhang =", opt$trimOverhang, "\n")
  }
  
  # Learn error model for forward reads.
  err_reverse <- learnErrors(fls = reverse_in, 
                             nreads = opt$num_learn, 
                             multithread = multithread_opt, 
                             randomize=opt$randomize)
  
  # Loop over all samples.
  for(sample in forward_samples) {
    
    if(opt$verbose) {
      cat("Processing: ", sample, "\n")
    }
    
    # Dereplicate reads.
    derep_forward <- derepFastq(forward_in[[sample]], n=opt$num_derep, verbose=opt$verbose)
    derep_reverse <- derepFastq(reverse_in[[sample]], n=opt$num_derep, verbose=opt$verbose)
    
    # Run inference and add inferred sequences to output list.
    dada_out_forward <- dada(derep_forward, err=err_forward, multithread=multithread_opt)
    dada_out_reverse <- dada(derep_reverse, err=err_reverse, multithread=multithread_opt)
    
    # Merge forward and reverse reads together
    seq_variants[[sample]] <- mergePairs(dadaF=dada_out_forward,
                                         derepF=derep_forward,
                                         dadaR=dada_out_reverse,
                                         derepR=derep_reverse,
                                         minOverlap=opt$minOverlap,
                                         maxMismatch=opt$maxMismatch,
                                         returnRejects=opt$returnRejects,
                                         justConcatenate=opt$justConcatenate,
                                         trimOverhang=opt$trimOverhang,
                                         verbose=opt$verbose)
    
    # Add read counts to log dataframe for debug purposes.
    log_counts[sample,] <- c(sample, sum(derep_forward$uniques), sum(derep_reverse$uniques), sum(dada_out_forward$denoised),
                            sum(dada_out_reverse$denoised), sum(seq_variants[[sample]]$abundance))
  }
}

if(opt$plot_errors) {
  
  if(opt$verbose){
    cat("Plotting estimated error models.\n\n")
  }
  
 if(opt$single){
   pdf("estimated_err.pdf", width=7, height=7)
   plotErrors(err_forward, nominalQ=TRUE)
   dev.off()
 } else{
   pdf("estimated_forward_err.pdf", width=7, height=7)
   plotErrors(err_forward, nominalQ=TRUE) 
   dev.off()
   
   pdf("estimated_reverse_err.pdf", width=7, height=7)
   plotErrors(err_reverse, nominalQ=TRUE) 
   dev.off()
 }
}

# Write out sequence table as RDS.
seqtab <- makeSequenceTable(seq_variants)
saveRDS(seqtab, opt$output)

# Write out logfile.
log_counts$tabled <- rowSums(seqtab)

write.table(x = log_counts, file = opt$log, quote = FALSE, sep="\t",
            col.names = TRUE, row.names = FALSE)
