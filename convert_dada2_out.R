#!/usr/bin/env Rscript

# Read in package to read in command-line options.
suppressMessages(library("optparse"))

version <- "1.0"

option_list <- list(
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input RDS file of dada2 sequence table (required)." , 
              metavar="path"),
  
  make_option(c("-b", "--biom"), type="character", default="seqtab.biom.tsv",
              help="Path to output BIOM file in TSV format (required).", 
              metavar="path"),
  
  make_option(c("-f", "--fasta"), type="character", default="seqtab.fasta",
              help="Path to output FASTA file (required).", 
              metavar="path"),
  
  make_option(c("--taxa_in"), type="character", default=NULL,
              help="Optional path to RDS containing table of taxa labels for dada2 variants", 
              metavar="path"),
  
  make_option(c("--taxa_out"), type="character", default="taxa_metadata.txt",
              help="Path for output of taxonomy metadata (default: taxa_metadata.txt).", 
              metavar="path"),
  
  make_option(c("--verbose"), action = "store_true", type="logical", default=FALSE,
              help="Print out more detailed information.", metavar = "boolean"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
  
)

opt_parser <- OptionParser(
                   option_list=option_list, 
                   usage = "%prog [options] -i seqtab_final.rds -b seqtab.biom -f seqtab.fasta",
                   description = paste(
                     "Script to convert dada2 sequence table to a BIOM (in legacy TSV) and FASTA file.\n",
                     "You can also specify an RDS file containing the taxa labels to create",
                     "a BIOM observation metadata table for taxonomy", sep="")
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
  stop("path to input sequence table RDS needs to be set.")
}

if(is.null(opt$biom)) {
  stop("path to output biom table needs to be set.")
}

if(is.null(opt$fasta)) {
  stop("path to output fasta file needs to be set.")
}

# Read in other required packages.
if(opt$verbose){
  library("ShortRead")
} else {
  suppressMessages(library("ShortRead"))
}

in_seqtab <- readRDS(opt$input)

seqs <- colnames(in_seqtab)
ids_study <- paste("seq", 1:ncol(in_seqtab), sep = "_")

colnames(in_seqtab) <- ids_study

if(opt$verbose){
  cat("\nWriting sequences to", opt$fasta, "\n\n")
}

# Write out fasta file.
writeFasta(ShortRead(sread = DNAStringSet(seqs), id = BStringSet(ids_study)), 
           file = opt$fasta)

if(opt$verbose){
  cat("Writing output biom table to", opt$biom, "\n\n")
}

# Transpose sequence table and add ids as new column.
in_seqtab_t <- as.data.frame(t(in_seqtab))
orig_col <- colnames(in_seqtab_t)
in_seqtab_t$ids <- rownames(in_seqtab_t)
in_seqtab_t <- in_seqtab_t[, c("ids", orig_col)]

col_df <- data.frame(matrix(c("#OTU ID", orig_col), nrow=1, ncol=ncol(in_seqtab_t)),
                     stringsAsFactors = FALSE)

colnames(col_df) <- colnames(in_seqtab_t)

# Write out BIOM file in legacy TSV format.
write.table(x = rbind(col_df, in_seqtab_t),
            file = opt$biom,
            quote = FALSE,
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE)

# If taxa file specified then also read that in.
if(!is.null(opt$taxa_in)){
  
  in_taxa <- readRDS(opt$taxa_in)
  
  # Check that dada2 variants are in same order as sequence table.
  if(!identical(rownames(in_taxa), seqs)){
    stop("sequences in taxa and sequence table are not ordered the same.")
  }
  
  rownames(in_taxa) <- ids_study

  # Replace all NA values with "Unclassified".
  in_taxa[is.na(in_taxa)] <- "Unclassified"
  
  taxa_combined <- apply(in_taxa, 1, function(x) paste(x, collapse=";"))
  
  taxa_out <- data.frame(names(taxa_combined), taxa_combined)
  colnames(taxa_out) <- c("#OTU ID", "taxonomy")
  
  if(opt$verbose){
    cat("Writing taxonomy observation metadata to", opt$taxa_out, "\n\n")
  }
  
  write.table(x = taxa_out, file = opt$taxa_out, quote = FALSE, sep="\t",
              col.names = TRUE, row.names = FALSE)
  
} else {
  
  if(opt$verbose){
    cat("No taxonomy file specified so not making taxonomy observation file.\n")
  }
   
}
