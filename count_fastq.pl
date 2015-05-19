#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use List::Util 'sum';

my ($parallel,$help);
my $res = GetOptions("parallel:i"=>\$parallel,
		     "help"=>\$help,
    )or pod2usage(2);

pod2usage(-verbose=>2) if $help;

my $cpu_count=0;
#if the option is set
if(defined($parallel)){
    #option is set but with no value then use the max number of proccessors
    if($parallel ==0){
	#load this module dynamically
	eval("use Sys::CPU;");
	$cpu_count=Sys::CPU::cpu_count();
    }else{
	$cpu_count=$parallel;
    }
}
my $pm = new Parallel::ForkManager($cpu_count);

my @files=@ARGV;


my $gzipped=0;
foreach my $input_file (@files){
    my $pid = $pm->start and next; 

    my ($file,$dir,$suffix)=fileparse($input_file, qr/\.[^.]*/);
    my $IN;
    if($suffix eq '.gz'){
	open($IN, "gunzip -c $input_file |") || die "can't open pipe to $input_file";
    }else{
	open($IN, '<',$input_file) || die "can't read $input_file";	
    }
    my $i=3;
    my $seq_count=0;
    my @lengths;
    while(<$IN>){
	#read every 4th line
	unless($i++ %4){
	    $seq_count++;
	    push @lengths,length($_);
	}
    }
    print join("\t",$file,$seq_count,sprintf('%d',sum(@lengths)/scalar(@lengths))),"\n";

    $pm->finish;
}

#Wait for all samples to be processed
$pm->wait_all_children;
