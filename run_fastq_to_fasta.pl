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

pod2usage($0.': You must provide a list of fastq files to be converted.') unless @files;

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
    my $out_file=$out_dir.'/'.$file_name.'.fasta';

    my $cmd="fastq_to_fasta -i $file -o $out_file";
    print $cmd,"\n";
    system($cmd);
    $pm->finish;
}

$pm->wait_all_children;

__END__

=head1 Name

run_fastq_to_fasta.pl - A wrapper of the FASTX toolkit command "fastq_to_fasta"

=head1 USAGE

run_fastq_to_fasta [-p [<# proc>] -o <out_dir> -h] <list of fastq files>

=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all FASTA output files.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_fastq_to_fasta.pl> This script simplifies running the FASTX toolkit "fastq_to_fasta" command by allowing multiple FASTQs to be converted with 1 command.

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

Minor updates by Gavin Douglas (gavin.douglas@dal.ca)

=cut

