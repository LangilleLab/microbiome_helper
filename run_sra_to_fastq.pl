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

pod2usage($0.': You must provide a list of sra files to be converted.') unless @files;

#make output directory 
system("mkdir -p $out_dir");

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

foreach my $file (@files){
    my $pid = $pm->start and next; 
    my ($file_name,$dir)=fileparse($file, qr/\.[^.]*/);
    my $out_file1=$out_dir.'/'.$file_name.'_1.fastq';
    my $out_file2=$out_dir.'/'.$file_name.'_2.fastq';

    unless(-e $out_file1.'.gz' && -e $out_file2.'.gz'){
	my $cmd="fastq-dump -F --split-files -O $out_dir $file";
	print $cmd,"\n";
	system($cmd);
	#note that the fastq-dump has a gzip output option but I found that this resulted in corupt output files.
	#therefore, just compress here instead
	my $gzip_cmd="gzip $out_file1 $out_file2";
	print $gzip_cmd,"\n";
	system($gzip_cmd);
    }   
    $pm->finish;
}

$pm->wait_all_children;

__END__

=head1 Name

run_sra_to_fastq.pl - A simple wrapper for the fastq-dump command for converting sra files to fastq

=head1 USAGE

run_sra_to_fastq.pl [-p [<# proc>] -o <out_dir> -h] <list of sra files>

run_sra_to_fastq.pl *.sra

#Specify alternate location for output files (instead of default current directory)

run_sra_to_fastq.pl -o fastq_files *.sra

#Run in parallel and use all CPUs

run_sra_to_fastq.pl *.sra -p

#Run in parallel limit to only 2 CPUs

run_sra_to_fastq.pl *.sra -p 2


=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all the output files.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_sra_to_fastq.pl> This script converts a directory of sra files into fastq files and gzips the output files as it proceeds


=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

