#!/usr/bin/perl
#Created by: Morgan Langille
#metaphlan_to_stamp.pl merged_metaphlan_output.txt > stamp_profile.spf

use warnings;
use strict;

my $header=<>;

chomp($header);

my @cols=split(/\t/,$header);

shift(@cols);

print join("\t","Kingdom","Phylum","Class","Order","Family","Genus","Species",@cols),"\n";

while(<>){
    chomp;
    my @data=split;
    my @taxonomy= split(/\|/,shift(@data));
    for my $x (0..6){
	unless (defined($taxonomy[$x])){
	    $taxonomy[$x]='Unclassified';
	}
    }
    print join("\t",@taxonomy,@data),"\n";
}
