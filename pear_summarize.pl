#!/usr/bin/perl

#pear_summarize.pl pear.log

use warnings;
use strict;
use List::Util qw(min max sum);

sub mean {
    return sum(@_)/@_;
}

my @samples;

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
	push (@samples, [$assembled_file,$assembled_percent,$discarded_percent,$unassembled_percent]);
    }
}



#sort in percent assembled in assending order
my @sorted_samples = sort {$a->[1] <=> $b->[1]} @samples;

unshift @sorted_samples,['Max',sprintf("%.3f",max(map{$_->[1]}@samples)),sprintf("%.3f",max(map{$_->[2]}@samples)),sprintf("%.3f",max(map{$_->[3]}@samples))];
unshift @sorted_samples,['Mean',sprintf("%.3f",mean(map{$_->[1]}@samples)),sprintf("%.3f",mean(map{$_->[2]}@samples)),sprintf("%.3f",mean(map{$_->[3]}@samples))];
unshift @sorted_samples,['Min',sprintf("%.3f",min(map{$_->[1]}@samples)),sprintf("%.3f",min(map{$_->[2]}@samples)),sprintf("%.3f",min(map{$_->[3]}@samples))];

#print header
print join("\t","ID","Assembled","Discarded","Unassembled"),"\n";

#print out all the data
foreach my $sample (@sorted_samples){
    print join("\t",@$sample),"\n";
}
