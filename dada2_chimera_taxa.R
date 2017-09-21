#!/usr/bin/env Rscript

# Read in package to read in command-line options.
library("optparse")

version <- "1.0"

option_list <- list(

  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Location of input .RDS file (required)." , metavar="path"),
  
  make_option(c("-r", "--refFasta"), type="character", default=NULL,
              help="Location of reference fasta to use for taxa assignment (required).", 
              metavar="path"),
  
  make_option(c("-s", "--ref_species"), type="character", default=NULL,
              help="Location of reference fasta to use for species assignment (required).", 
              metavar="path"),
  
  make_option(c("--seed"), type="integer", default=NULL,
              help="Random seed to make command reproducible (required).", metavar = "integer"),
  
  make_option(c("-t", "--threads"), type="integer", default=1,
              help="Number of threads to use (default: 1).", metavar="integer"),
  
  make_option(c("--skip_species"), action = "store_true", type="logical", default=FALSE,
              help=paste("Flag to indicate that the species assignment step should be skipped.",
                         "This step can often be the slowest step. (default: FALSE).", sep=" "),
              metavar = "boolean"),
  
  make_option(c("--count_out"), type="character", default="seqtab_final.rds",
              help=paste("Output RDS file for sequence table after chimera filtering",
                         "(default: \"seqtab.rds\").", sep=" "), metavar="path"),
  
  make_option(c("--tax_out"), type="character", default="tax_final.rds",
              help=paste("Output RDS file for sequence table after chimera filtering",
                         "(default: \"tax_final.rds\").", sep=" "), metavar="path"),
  
  make_option(c("--chimera_method"), type="character", default="consensus",
              help=paste("Method for chimera checking, one of \"consensus\", \"pooled\",",
                         "or \"per-sample\". See dada2 manual for details.",
                         "(default: \"consensus\").", sep=" ")),
  
  make_option(c("--minFoldParentOverAbundance"), type="numeric", default=1,
              help=paste("Sequences with this fold-abundance greater than a given sequence", 
                         "can be considered as \"parents\" (default: 1)", sep=" "), 
              metavar="numeric"),
  
  make_option(c("--minParentAbundance"), type="numeric", default=8,
              help=paste("Sequences with this abundance", 
                         "can be considered as \"parents\" (default: 8)", sep=" "), 
              metavar="numeric"),
  
  make_option(c("--allowOneOff"), action = "store_true", type="logical", default=TRUE,
              help=paste("Flag to indicate that sequences with one mismatch or indel",
                         "to an exact bimera are also taken to be bimeras (default: TRUE).", sep=" "),
              metavar = "boolean"),
  
  make_option(c("--minOneOffParentDistance"), type="numeric", default=4,
              help=paste("When flagging one-off bimeras only sequences with at least this many", 
                         "many mismatches to potential bimera are considered \"parents\" (default: 4)", 
                         sep=" ")), 
              
 make_option(c("--maxShift"), type="numeric", default=16,
             help=paste("Max shift allowed when aligning sequences to potential \"parents\"", 
                        "(default: 16)", sep=" ")), 

 make_option(c("--minSampleFraction"), type="numeric", default=0.9,
             help=paste("Fraction of samples in which a sequence must be flagged as a bimera", 
                        "to be classified as a bimera (default: 0.9)", sep=" ")),               
 
 make_option(c("--ignoreNNegatives"), type="integer", default=1,
             help=paste("Number of unflagged samples to ignore when evaluating whether", 
                        "sequence was flagged in bimera in more than minSampleFraction",
                        "samples (default: 1)", sep=" ")), 
 
 make_option(c("--minBoot"), type="numeric", default=50,
             help=paste("Min bootstrap confidence for assigning taxa label", 
                        "(default: 50)", sep=" ")), 
 
 make_option(c("--tryRC"), action = "store_true", type="logical", default=FALSE,
             help=paste("Flag to indicate that the reverse-complement of sequences",
                        "should be used for taxa assignment if it is a better match",
                        "to the reference database (default: FALSE).", sep=" "),
             metavar = "boolean"),
 
 make_option(c("--allowMultiple"), action = "store_true", type="logical", default=FALSE,
             help=paste("Flag to indicate that if multiple different species are exact",
                        "matches then these should all be returned (default: FALSE).", sep=" "),
             metavar = "boolean"),
 
  make_option(c("--log"), type="character", default="dada2_nonchimera_counts.txt",
              help="Output logfile for read count table (default: dada2_nonchimera_counts.txt).", 
              metavar="path"),
  
  make_option(c("--verbose"), action = "store_true", type="logical", default=FALSE,
              help="Write out status messages (default: FALSE).", metavar = "boolean"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")

)

opt_parser <- OptionParser(option_list=option_list, 
                           usage = "%prog [options] -i seqtab.rds -r /path/to/ref.fa -s /path/to/species_ref.fa --seed 123",
                           
                           description = paste("\nThis is a wrapper for the DADA2 chimera and tax assignment step that is", 
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

# Check if path to input RDS set, if not stop job.
if(is.null(opt$input)) {
  stop("path to input RDS needs to be set.")
}

# Check that refFasta set.
if(is.null(opt$refFasta)) {
  stop("Path to refFasta needs to be set.")
}

# Check that ref_species set if doing species assignment.
if(!opt$skip_species & is.null(opt$refFasta)) {
  stop("Path to refFasta needs to be set.")
}

# Check that random seed set.
if(is.null(opt$seed)) {
  stop("random seed needs to be set.")
}

# Check that chimera method is one of 3 possible choices.
if(! opt$chimera_method %in% c("consensus", "pooled", "per-sample")) {
 stop(paste("--chimera_method needs to be one of \"consensus\",",
            "\"pooled\", or \"per-sample\"", sep=" "))
}

# Set multithread option (FALSE if no, otherwise give # core).
if (opt$threads > 1) {
  multithread_opt <- opt$threads
} else {
  multithread_opt <- FALSE
}

# Load in input RDS file.
in_seqtab <- readRDS(opt$input)

# Read in and print out DADA2 version.
library("dada2")

if(opt$verbose) {
        cat("Running DADA2 version:", as.character(packageVersion("dada2")), "\n\n",
        "Running removeBimeraDenovo with below options.\n",
        "method =", opt$chimera_method, "\n",
        "minFoldParentOverAbundance =", opt$minFoldParentOverAbundance, "\n",
        "minParentAbundance =", opt$minParentAbundance, "\n",
        "allowOneOff =", opt$allowOneOff, "\n",
        "minOneOffParentDistance =", opt$minOneOffParentDistance, "\n",
        "maxShift =", opt$maxShift, "\n",
        "minSampleFraction =", opt$minSampleFraction, "\n",
        "ignoreNNegatives =", opt$ignoreNNegatives, "\n\n")
}

seqtab_nochim <- removeBimeraDenovo(in_seqtab, 
                                    method=opt$chimera_method, 
                                    minFoldParentOverAbundance=opt$minFoldParentOverAbundance,
                                    minParentAbundance=opt$minParentAbundance,
                                    allowOneOff=opt$allowOneOff,
                                    minOneOffParentDistance=opt$minOneOffParentDistance,
                                    maxShift=opt$maxShift,
                                    minSampleFraction=opt$minSampleFraction,
                                    ignoreNNegatives=opt$ignoreNNegatives,
                                    multithread=multithread_opt,
                                    verbose=opt$verbose)

if(opt$verbose) {
  cat("Running assignTaxonomy with below options.\n",
      "Random seed:", as.character(opt$seed), "\n\n",
      "refFasta =", opt$refFasta, "\n",
      "minBoot =", opt$minBoot, "\n",
      "tryRC =", opt$tryRC, "\n\n")
}

# Set random seed for reproducibility.
set.seed(opt$seed)

# Run taxonomy assignment.
taxa <- assignTaxonomy(seqs=seqtab_nochim, 
                       refFasta=opt$refFasta,
                       minBoot=opt$minBoot,
                       tryRC=opt$tryRC,
                       multithread=multithread_opt,
                       verbose=opt$verbose)

if(! opt$skip_species) {

  if(opt$verbose) {
    cat("\n\nRunning addSpecies with below options.\n",
        "refFasta =", opt$ref_species, "\n",
        "allowMultiple =", opt$allowMultiple, "\n\n")
  }
  
  # Add exact species assignment where possible.
  taxa_species <- addSpecies(taxtab=taxa, 
                             refFasta=opt$ref_species, 
                             allowMultiple = opt$allowMultiple, 
                             verbose = opt$verbose)
  
  # Write out taxa + species table.
  saveRDS(taxa_species, opt$tax_out)
  
} else {
  
  if(opt$verbose) {
    cat("\n\nSkipping addSpecies step.\n\n")
  }
 
  # Write out the taxa count table without species.
  saveRDS(taxa, opt$tax_out)
  
}

# Write out count table as RDS.
saveRDS(seqtab_nochim, opt$count_out)


# Track read and variants counts before and after chimera checking.

log_counts <- as.data.frame(
               cbind(rowSums(in_seqtab),
               rowSums(seqtab_nochim),
               rowSums(in_seqtab != 0),
               rowSums(seqtab_nochim != 0)))

saveRDS(in_seqtab,file = "tmp1.rds")
saveRDS(seqtab_nochim,file = "tmp2.rds")
saveRDS(in_seqtab,file = "tmp3.rds")
saveRDS(seqtab_nochim,file = "tmp4.rds")

colnames(log_counts) <- c("input_reads", "nonchimera_reads",
                          "input_variants", "nonchimera_variants")

log_counts$sample <- rownames(in_seqtab)

log_counts <- log_counts[,c("sample", "input_reads", "nonchimera_reads",
                            "input_variants", "nonchimera_variants")]

write.table(x = log_counts, file = opt$log, quote = FALSE, sep="\t",
            col.names = TRUE, row.names = FALSE)
