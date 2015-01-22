#!/usr/bin/perl

#humann_to_stamp.py 04b-hit-keg-mpm-cop-nul-nve-nve.txt > hummann_modules.spf
#humann_to_stamp.py 04b-hit-keg-mpt-cop-nul-nve-nve.txt > hummann_pathways.spf
#humann_to_stamp.py 01b-hit-keg-cat.txt > hummann_kos.spf

use warnings;
use strict;

$_=<>;
chomp;
my @header=split(/\t/,$_);
shift @header;
foreach my $x(0..$#header){
    $header[$x]=~ s/-hit-keg-(.*)//; 
}
print join("\t",@header),"\n";

my $x=0;
while(<>){
    if($x==0){
	if(/^Richness/){
	    $x=1;
	    next;
	}
    }else{
	chomp;
	my @values=split(/\t/,$_);
	shift @values;
	print join("\t",@values),"\n";
    }
    
}
