#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my ($parallel,$help);
my $out_dir='./';
my $index_dir = "/home/shared/bowtiedb/hg19";
my $res = GetOptions("out_dir=s" => \$out_dir,
		     "db=s" => \$index_dir,
		     "parallel:i"=>\$parallel,
		     "help"=>\$help,
    )or pod2usage(2);

pod2usage(-verbose=>2) if $help;

my @files=@ARGV;

pod2usage($0.': You must provide a list of fasta files to be screened for human sequences.') unless @files;

#make output directory 
system("mkdir -p $out_dir");

my $cpu_count=1;
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

my %paired_files;
foreach my $file (@files){
    my ($file_name,$dir,$suffix)=fileparse($file, qr/\.[^.]*/);

    my $out_file=$out_dir.'/'.$file_name."_screened".$suffix;

    my $cmd;
    if($suffix eq '.fastq'){
	#Can add --local option to make this more sensitive
	$cmd="bowtie2 -x $index_dir -p $cpu_count -U $file --un $out_file >/dev/null";
    }else{
	$cmd="bowtie2 -x $index_dir -p $cpu_count -f $file --un $out_file >/dev/null";
    }
    print $cmd,"\n";
    system($cmd);
}

__END__

=head1 Name

run_human_filter.pl - A simple wrapper for human_filter to screen for human sequences in metagenomic data

=head1 USAGE

run_human_filter.pl [-p [<# proc>] -o <out_dir> -d <path to index files> -h] <list of fastq or fasta files>

E.g.

run_human_filter.pl sample1_assembled.fastq sample2_assembled.fastq

#Shorter way to do the same thing

run_human_filter.pl *.fastq 

#Same as above, but specify a different bowtie index database

run_human_filter.pl *.fastq -d /your_path/prefix

#Specify alternate location for output files (instead of default current directory)

run_human_filter.pl -o screened_reads *.fastq

#Run in parallel and use all CPUs

run_human_filter.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

run_human_filter.pl *.fastq -p 2




=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all output files.

=item B<-d, --db <path>>

Path to human genome bowtie index files. By default is "/home/shared/bowtiedb/hg19".
Note that the last part of the path is actually the prefix of the index files 
(e.g. the files in /home/shared/bowtiedb/ are: hg19.1.bt2  hg19.2.bt2  hg19.3.bt2  hg19.4.bt2  hg19.rev.1.bt2  hg19.rev.2.bt2). 

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_human_filter.pl> This script allows for more automated running of the human_filter program on multiple fastq files. 
 
Before use make sure you have installed the "human_filter.pl" (along with corresponding databases) program so it is accesible from your PATH

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

Updated 30/03/2016 by Gavin Douglas  E<lt>gavin.douglas@dal.caE<gt>

=cut

