#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my ($parallel,$help);
my $out_dir='./';
my $res = GetOptions("out_dir=s" => \$out_dir,
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
my $pm = new Parallel::ForkManager($cpu_count);

my %paired_files;
foreach my $file (@files){
    my $pid = $pm->start and next; 
    my $cmd="deconseq.pl -f $file -dbs hsref -i 98 -c 90 -out_dir $out_dir";
    print $cmd,"\n";
    system($cmd);
    $pm->finish;
}

#Wait for all samples to be processed
$pm->wait_all_children;

__END__

=head1 Name

run_deconseq.pl - A simple wrapper for deconseq to screen for human sequences in metagenomic data

=head1 USAGE

run_deconseq.pl [-p [<# proc>] -o <out_dir> -h] <list of fastq or fasta files>

E.g.

run_deconseq.pl sample1_assembled.fastq sample2_assembled.fastq

#Shorter way to do the same thing

run_deconseq.pl *.fastq

#Specify alternate location for output files (instead of default current directory)

run_deconseq.pl -o screened_reads *.fastq

#Run in parallel and use all CPUs

run_deconseq.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

run_deconseq.pl *.fastq -p 2


=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all output files.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_deconseq.pl> This script allows for more automated running of the deconseq program on multiple fastq files. 
 
Before use make sure you have installed the "deconseq.pl" (along with corresponding databases) program so it is accesible from your PATH

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

