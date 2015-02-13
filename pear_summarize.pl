#!/usr/bin/perl

#pear_summarize.pl pear.log

use warnings;
use strict;

#print header
print join("\t","Assembled","Discarded","Unassembled","File"),"\n";

while(<>){
    chomp;
    if(/Assembled reads/){
	my $assembled_string=$_;
	my $discarded_string=<>;
	my $unassembled_string=<>;
	my $assembled_file_string=<>;
	my ($assembled_percent) = $assembled_string =~ /(\d+\.\d+)\%/;
	my ($discarded_percent) = $discarded_string =~ /(\d+\.\d+)\%/;
	my ($unassembled_percent) = $unassembled_string =~ /(\d+\.\d+)\%/;
	my ($assembled_file) = $assembled_file_string =~ /(\w+)\.assembled/;
	print join("\t",$assembled_percent,$discarded_percent,$unassembled_percent,$assembled_file),"\n";
    }
}
