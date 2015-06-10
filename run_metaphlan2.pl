#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my $metaphlan_dir='/usr/local/metaphlan2/';

#location to store intermediate metaphlan files
my $metaphlan_out_dir='./metaphlan_out/';

my $bowtie='sensitive-local';
my $align_len=50;

my ($final_out_file,$parallel,$help);
my $res = GetOptions("output=s" => \$final_out_file,
		     "location=s"=> \$metaphlan_dir,
		     "align_len=i"=>\$align_len,
		     "parallel:i"=>\$parallel,
		     "bowtie:s"=>\$bowtie,
		     "help"=>\$help,
    )or pod2usage(2);

pod2usage(-verbose=>2) if $help;

pod2usage($0.': You must specify an output file.') unless defined $final_out_file;

my $metaphlan_script=$metaphlan_dir.'metaphlan2.py';
my $metaphlan_db=$metaphlan_dir.'db_v20/mpa_v20_m200';
my $metaphlan_pkl=$metaphlan_dir.'db_v20/mpa_v20_m200.pkl';
my $metaphlan_merge=$metaphlan_dir.'utils/merge_metaphlan_tables.py';

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
system("mkdir -p $metaphlan_out_dir");

my @files=@ARGV;

my $format='multifasta';

my %paired_files;
my $gzipped=0;
foreach(@files){
    my ($file,$dir,$suffix)=fileparse($_, qr/\.[^.]*/);
    if($suffix eq '.gz'){
	$gzipped=1;
    }
    if($suffix =~ /fastq/ || $file =~ /fastq/){
	$format='multifastq';
    }

    my $name=$file;
    if($file =~ /(.+)_R[1|2]_/){
	$name=$1;
    }
    push(@{$paired_files{$name}},$_);
}

my @out_files = map ($metaphlan_out_dir.$_, keys %paired_files);

foreach my $name (keys %paired_files){
    my $pid = $pm->start and next; 
    my $cat;
    if ($gzipped){
	$cat='zcat';
    }else{
	$cat='cat';
    }
    my $out_file=$metaphlan_out_dir.$name;
    my $cmd=join(' ',$cat,@{$paired_files{$name}});
    $cmd.=" | $metaphlan_script  --input_type $format --mpa_pkl $metaphlan_pkl --bt2_ps $bowtie --min_alignment_len $align_len --bowtie2db $metaphlan_db --no_map > $out_file";
    print $cmd,"\n";
    system($cmd);
    $pm->finish;
}

#Wait for all samples to be processed
$pm->wait_all_children;

#merge metaphlan output
my $merge_cmd=$metaphlan_merge.' '.join(' ',@out_files).' > '.$final_out_file;
print $merge_cmd,"\n";
system($merge_cmd);

__END__

=head1 Name

run_metaphlan2.pl - Provide a simpler way to run metaphlan

=head1 USAGE

run_metaphlan2.pl [-b <bowtie setting> -a <min_align_len> -l <metaphlan_dir> -p [<# proc>] -h] -o out.txt <list of fastq files>

E.g.

run_metaphlan2.pl -o out.txt sample1.fastq sample2.fastq sample3.fastq

#shorter way of writing the same thing

run_metaphlan2.pl -o out.txt *.fastq

#Run in parallel and use all CPUs

run_metaphlan2.pl -o out.txt *.fastq -p

#Run in parallel limit to only 2 CPUs

run_metaphlan2.pl -o out.txt *.fastq -p 2

#fastq files can be gzipped (note: all files must be either gzipped or not. Can't be a mix)

run_metaphlan2.pl -o out.txt *.fastq.gz

#input files can be in fasta format in addition to fastq (optionally gzipped as well)

run_metaphlan2.pl -o out.txt *.fasta

#paired end files can be handled by concatentating them on the fly (files must have "_R1_" and "_R2_" within the file name)

run_metaphlan2.pl -o out.txt sample1_R1_001.fastq.gz sample1_R2_001.fastq.gz sample2_R1_001.fastq.gz sample2_R2_001.fastq.gz


=head1 OPTIONS

=over 4

=item B<-o, --output <file>>

Mandatory. The name of the file for the merged data to be written to.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-a, --align_len <min_align_len>>

This sets the minimum alignment length for sequences to match between the reads and the marker database.

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_metaphlan2.pl> This script allows for easier running of the metaphlan pipeline for taxonomy assingment on metagenomic data. In particular it automates the running of multiple metagenomic samples and merges the data into a single output table. It handles combining paired end data by contatentating the files on the fly. It also handles gzipped fastq files without creating uncompressed intermediate files. 

Before use make sure you have installed Bowtie2 and the metaphlan package.

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

