#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use List::Util 'sum';

my ($parallel,$len_cutoff,$help);
my $out_dir='./';
my $res = GetOptions("out_dir=s" => \$out_dir,
		     "parallel:i"=>\$parallel,
		     "len_cutoff:i"=>\$len_cutoff,
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

#make output directory 
system("mkdir -p $out_dir");

my $gzipped=0;
foreach my $input_file (@files){
    my $pid = $pm->start and next; 

    my ($file,$dir,$suffix)=fileparse($input_file, qr/\.[^.]*/);

    print "Filtering: $input_file\n";
    my $out_file=$out_dir.'/'.$file.".filtered".$suffix;

    open (my $OUT,'>',$out_file) || die "Can't open file for writing: $out_file";

    my $IN;
    if($suffix eq '.gz'){
	open($IN, "gunzip -c $input_file |") || die "can't open pipe to $input_file";
    }else{
	open($IN, '<',$input_file) || die "can't read $input_file";	
    }
    my $i=3;
    my @lengths;
    while(my $id=<$IN>){
	my $seq=<$IN>;
	my $strand=<$IN>;
	my $quality=<$IN>;
	if(!($seq =~ /n|N/)){
	    #subtract one since the sequence always has a line return character at the end.
	    if(!$len_cutoff || length($seq)-1 >= $len_cutoff){
		print $OUT $id;
		print $OUT $seq;
		print $OUT $strand;
		print $OUT $quality;
	    }
	}
	    
    }
    $pm->finish;
}

#Wait for all samples to be processed
$pm->wait_all_children;
