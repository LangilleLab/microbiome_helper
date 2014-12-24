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

pod2usage($0.': You must provide a list of fastq files to be merged.') unless @files;

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
    my ($file_name,$dir)=fileparse($file);
    if($file_name =~ /(.+)_R([1|2])_/){
	$paired_files{$1}[$2-1]=$file;
    }else{
	warn "Input file does not contain '_R1_' or '_R2_' in name: $file";
    }
}

foreach my $name (keys %paired_files){
    unless(defined($paired_files{$name}[0]) && defined($paired_files{$name}[1])){
	warn "Couldn't find matching paired end files for file starting with: $name";
	next;
    }
    my $out_file=$out_dir.'/'.$name;
    my $cmd="pear -f $paired_files{$name}[0] -r $paired_files{$name}[1] -j $cpu_count -o $out_file";
    print $cmd,"\n";
    system($cmd);
}

__END__

=head1 Name

run_pear.pl - A simple wrapper for PEAR to stich paired-end reads

=head1 USAGE

run_pear.pl [-p [<# proc>] -o <out_dir> -h] <list of fastq files>

E.g.

#Note: files must have "_R1_" and "_R2_" within the file name

run_pear.pl sample1_R1_001.fastq sample1_R2_001.fastq sample2_R1_001.fastq sample2_R2_001.fastq

#Shorter way to do the same thing

run_pear.pl *.fastq

#Specify alternate location for output files (instead of default current directory)

run_pear.pl -o stitched_reads *.fastq

#Run in parallel and use all CPUs

run_pear.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

run_pear.pl *.fastq -p 2


=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all PEAR output files.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_pear.pl> This script allows for more automated running of the PEAR program on multiple fastq files. PEAR is used to stitch (or assemble) paired end reads together. The assumption is made that the paired end files have the same name with the forward reads being indicated by "_R1_" and the reverse being "_R2_". The script also allows the use of multiple threads.
 
Before use make sure you have installed the "pear" program so it is accesible from your PATH

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

