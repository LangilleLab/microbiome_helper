#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;


my ($gzip,$write_only,$parallel,$help);
my $out_dir='./';
my $res = GetOptions("out_dir=s" => \$out_dir,
		     "parallel:i"=>\$parallel,
		     "gzip"=>\$gzip,
		     "write_only"=>\$write_only,
		     "help"=>\$help,
    )or pod2usage(2);

pod2usage(-verbose=>2) if $help;

my @files=@ARGV;

pod2usage($0.': You must provide a list of files to be concatenated.') unless @files;


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

#create output directory
system("mkdir -p $out_dir");

my %paired_files;
foreach my $file (@files){
    my ($file_name,$dir,$suffix)=fileparse($file,qr/(\.fasta|\.fastq)\.?[^.]*/);
    if($file_name =~ /(.+)_L00([\d])(_R[1|2])_/){
	my $name=$1;
	my $lane=$2;
	my $paired_end=$3;
	my $new_name=$name.$paired_end.$suffix;
	$paired_files{$new_name}[$lane-1]=$file;
    }else{
	warn "Input file does not contain '_R1_'and '_R2_' in name: $file";
    }
}


foreach my $name (keys %paired_files){
    my $pid = $pm->start and next; 

    my $cmd=join(' ','cat',@{$paired_files{$name}});
    
    my $out_file=$out_dir.'/'.$name;
    $cmd.= " > $out_file";
    print $cmd,"\n";
    unless($write_only){
	system($cmd);
    }
    $pm->finish;
}

__END__

=head1 Name

concat_lanes.pl - Concatenate multiple lanes of sequence files (fasta or fastq) (gzipped or not)

=head1 USAGE

concat_lanes.pl [-p [<# proc>] -h] -o <out_dir> <list of paired end files>

E.g.

#Note: Files must be named "L001_R1_", "L002_R1_", "L003_R1_", etc. 
#Results in lanes merged for both R1 and R2 separately

concat_lanes.pl sample1_L001_R1_001.fastq sample1_L002_R1_001.fastq sample2_L001_R2_001.fastq sample2_L002_R2_001.fastq

#Shorter way to do the same thing

concat_lanes.pl *.fastq

#Files can be gzipped

concat_lanes.pl *.fastq.gz

#Specify alternate location for output files (instead of default current directory)

concat_lanes.pl -o concatenated_reads *.fastq

#Run in parallel and use all CPUs

concat_lanes.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

concat_lanes.pl *.fastq -p 2

=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all  output files.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-w, --write_only>

Write the commands that will be run without actually running them.

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<concat_lanes.pl> 

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

