#!/usr/bin/perl
#Created by: Morgan Langille

#Note: metaphlan file must contain merged files
#Note2: 'all' taxonomy must be output in metaphlan

#metaphlan_to_stamp.pl merged_metaphlan_output.txt > stamp_profile.spf

use warnings;
use strict;

my $header=<>;

chomp($header);

my @cols=split(/\t/,$header);

shift(@cols);

print join("\t","Kingdom","Phylum","Class","Order","Family","Genus","Species",@cols),"\n";

LINE:
while(<>){
    chomp;
    my @data=split;
    my @taxonomy= split(/\|/,shift(@data));
    my $unclassified_flag=0;
    for my $x (0..6){
	if(defined($taxonomy[$x])){
	    if($taxonomy[$x] =~ /unclassified/){
		$taxonomy[$x]='unclassified';
		$unclassified_flag=1;
	    }
	}elsif($unclassified_flag ==1){
	    $taxonomy[$x]='unclassified';
	}else{
	    #not defined and not an unclassified parent so do not include in output
	    next LINE;
	}
    }
    print join("\t",@taxonomy,@data),"\n";
}

