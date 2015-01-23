#!/usr/bin/perl

#Creates a qiime map file based on a list of input filenames that can then be used with the script "add_qiime_labels.py"

#create_qiime_map.pl fasta_files/* >map.txt

use warnings;
use strict;
use File::Basename;

my @files=@ARGV;

#output header
print join("\t",'#SampleID','BarcodeSequence','LinkerPrimerSequence','FileInput'),"\n";

foreach my $file (@files){
    my ($file_name,$dir)=fileparse($file);
    my $sample_id=$file_name;

    #This attempts to remove sequencer labels that probably should belong in the sample_id
    if($file_name =~ /(_S\d+_L\d+.*)/){
	$sample_id =~ s/$1//;
    }
    #this removes all non alphanumeric characters and underscores and replaces them with a '.'
    $sample_id=~ s/[\W_]/\./g;
    print join("\t",$sample_id, '','',$file_name),"\n";
}


