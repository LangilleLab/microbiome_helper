#!/usr/bin/perl

use warnings;
use strict;

use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my ($parallel,$help);
my $out_dir = "screened_reads";
my $index_dir;
my $log = "screened_reads.log";
my $keep;

my $param_file;

my $res = GetOptions("out_dir=s" => \$out_dir,
		     "db=s" => \$index_dir,
		     "parallel:i"=>\$parallel,
		     "help"=>\$help,
		     "log=s"=>\$log,
		     "config_file=s"=>\$param_file,
		     "keep"=>\$keep,
    )or pod2usage(2);

pod2usage(-verbose=>2) if $help;

my @files=@ARGV;

my @log_files = ();

pod2usage($0.': You must provide a list of fasta files to be screened for contaminant sequences.') unless @files;

if ( ! $index_dir )	{	die "you need to specific a bowtie2 index directory with --db\n";	}

#check that all files exist:
foreach my $f ( @files )	{	
	if ( ! -e $f )	{	die "Stopping: file $f not found\n";	}
}

#if user specified parameter file
my $param = "";
if ( $param_file )	{
	open( 'CONFIG' , '<' , $param_file ) or die "cant open user specified parameter file $param_file\n";
	while( <CONFIG> )	{
		my @s = split( '\s+' , $_ );
		if ( ! exists $s[0] )	{	next	}	#blank line
		my $line = join( " " , @s );
		$param = $param . " " . $line;
	} close( 'CONFIG' );
}

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

foreach my $file (@files)	{

    my ($file_name,$dir,$suffix)=fileparse($file, qr/\.[^.]*/);

    my $out_file=$out_dir.'/'.$file_name."_screened".$suffix;

    my $cmd;

    my $log_out = $out_file;
    $log_out =~ s/$suffix$/.log/;

    push( @log_files , $log_out );

    if( $suffix eq '.fastq' )	{

	$cmd="bowtie2 -x $index_dir $param -p $cpu_count -U $file --un $out_file >/dev/null 2> $log_out";

    } else {

	$cmd="bowtie2 -x $index_dir $param -p $cpu_count -f $file --un $out_file >/dev/null 2> $log_out";

    }

    print $cmd,"\n";
    system($cmd);
}


open( 'LOG' , '>' , $log ) or die "cant create log summary $log\n";

print LOG "sample	input	aligned	percent_removed\n";

foreach my $l ( @log_files )	{

	my $total = 0;
	my $aligned = 0;
	my $percent = 0;

	open( 'TMP' , '<' , $l ) or die "cant open individual log file $l\n";
	while( <TMP> )	{
		
		if ( $_ =~ m/(\d+) reads; of these:/ )	{
			$total = $1;
		} elsif ( $_ =~ m/(\d+) .+ aligned exactly 1 time/ ) {
			$aligned += $1;
		} elsif ( $_ =~ m/(\d+) .+ aligned >1 times/ )  {
			$aligned += $1
		} elsif ( $_ =~ m/(\S+) overall alignment rate/ )	{
			$percent = $1;
		} else {}		

	} close( 'TMP' );
	
	my $sample = $l;
	$sample =~ s/\.log$//;

	print LOG "$sample	$total	$aligned	$percent\n";

	if ( ! $keep )	{	system( "rm $l" )	}

}

close( 'LOG' );

__END__

=head1 Name

run_contaminant_filter.pl - A simple wrapper for bowtie2 to screen out contaminant sequences in metagenomic data

=head1 USAGE

run_contaminant_filter.pl [-p [<# proc>] -o <out_dir> -d <path to index files> -l <log_file> -c <bowtie2 param file> -h --keep] <list of fastq or fasta files>

E.g.

run_contaminant_filter.pl sample1_assembled.fastq sample2_assembled.fastq

#Shorter way to do the same thing

run_contaminant_filter.pl *.fastq 

#Same as above, but specify a different bowtie index database

run_contaminant_filter.pl *.fastq -d /your_path/prefix

#Specify alternate location for output files (instead of default current directory)

run_contaminant_filter.pl -o screened_reads *.fastq

#Run in parallel and use all CPUs

run_contaminant_filter.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

run_contaminant_filter.pl *.fastq -p 2


=head1 OPTIONS

=over 4

=item B<-d, --db <path>>

Path to bowtie index files. Note that the last part of the path is actually the prefix of the index files 
(e.g. the files in /home/shared/bowtiedb/ are: hg19.1.bt2  hg19.2.bt2  hg19.3.bt2  hg19.4.bt2  hg19.rev.1.bt2  hg19.rev.2.bt2) and this would be specified with --db /home/shared/bowtiedb/hg19. 

=item B<-o, --out_dir <file>>
	
The name of the output directory to place all output files ("screened_reads" by default).
	
=item B<-l, --log <file>> 
	
The name of the summary logfile (default: "screened_reads.log"). 
	
=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-h, --help>

Displays the entire help documentation.

=item B<-c, --config_file>

Optional file that contains parameters to pass to bowtie2. There should be 1 parameter per line.

For example:

--local
--seed 100
-k 1

=item B<-k, --keep>

Flag that indicates log files for each individual command should not be deleted (helpful for troubleshooting).

=back

=head1 DESCRIPTION

B<run_contaminant_filter.pl> This is a wrapper script for running bowtie2 to screen out contaminant sequences (e.g. human sequences).
 
=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

Last updated 14 Dec 2016 by Gavin Douglas  E<lt>gavin.douglas@dal.caE<gt>

=cut

