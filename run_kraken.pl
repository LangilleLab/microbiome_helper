#!/usr/bin/perl

use warnings;
use strict;

use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $help; 
my $version = "0.1";
my $version_marker;

my $parallel = 1;
my $db ;
my $out_dir = "kraken_output";
my $log = "kraken_log.txt";
my $report = "kraken_mpa_report.txt";
my $keep;  

### delimiter used to identify sample name (either "_" or ".")
my $delimiter = "_";

### kraken parameters:
my $preload;
my $quick;
my $min_hits;
my $paired;
my $check_names;

### initialized variables not user-specified
my %s2f = (); ### samples as keys and fastqs as values (2 for PE and 1 for SE)
my @log_files = (); ### array to contain all log files, which will be looped over at the end
my @out_files = (); ### out files needed for kraken-mpa-report
my $kraken_options = ""; ### options to be passed to command-line

my $res = GetOptions("database|db=s"=>\$db,
		     "out_dir|o=s" => \$out_dir,
		     "log=s"=>\$log,
		     "thread:i"=>\$parallel,
		     "preload"=>\$preload,
		     "quick"=>\$quick,
		     "min-hits:i"=>\$min_hits,
		     "paired"=>\$paired,
		     "check-names"=>\$check_names,
		     "keep"=>\$keep,
		     "report=s"=>\$report,
		     "delimiter=s"=>\$delimiter,
		     "help|h"=>\$help,
		     "version|v"=>\$version_marker,
	  )	or pod2usage(2);

pod2usage(-verbose=>2) if $help;

if ( $version_marker )	{	print "version $version\n";	exit	}

### the "2>&1" syntax means output STDERR to same output as STDOUT (i.e. same output as file descripter 1).
if ( index( `kraken 2>&1` , "Need to specify input filenames!" ) == -1 )      { die "Stopping job, because \"kraken\" is not in your path.\n";	} 	

if ( ! defined $db )	{	die "need database directory\n";	}

if ( ( $delimiter ne "_" ) and ( $delimiter ne "." ) )	{	die "delimiter $delimiter needs to be \"_\" or \".\"\n";	}

my @files = @ARGV;

### link sample names to fastqs
&fastqs2sample( @files );

pod2usage($0.': You must provide a list of fastq files.') unless @files;

# check that all files exist: 
foreach my $f ( @files )	{	 
	if ( ! -e $f )	{ die "file $f doesn't exist\n";	}	
}

system( "mkdir -p $out_dir");

if (  $preload )	{	$kraken_options = $kraken_options . " --preload";	}
if (  $quick )	{	$kraken_options = $kraken_options . " --quick";	}
if (  defined $min_hits )	{	$kraken_options = $kraken_options . " --min-hits $min_hits";	}
if ( $paired )	{	$kraken_options = $kraken_options . " --paired";	}
if ( $check_names )	{	$kraken_options = $kraken_options . " --check-names";	}

foreach my $sample ( keys %s2f )	{
	
	my @fastqs = @{ $s2f{$sample} };

	my $out = "$out_dir" . "/" . $sample . "_out.txt";
	my $log = "$out_dir" . "/" . $sample . "_out.log";
			
	push( @out_files , $out );
	push( @log_files , $log );
			
	print STDERR "\nRunning:\nkraken --db $db $kraken_options --output $out @fastqs 2> $log\n";
	system( "kraken --db $db $kraken_options --output $out @fastqs 2> $log" );
	
}

### now combine all kraken output files into mpa report:
print STDERR "\nRunning:\nkraken-mpa-report --db $db @out_files > $report\n";
system( "kraken-mpa-report --header-line --db $db @out_files > $report" );

my $report_spf = $report; 
$report_spf =~ s/\.txt$/.spf/;

print STDERR "\nRunning:\nmetaphlan_to_stamp.pl $report > $report_spf\n";
system( "metaphlan_to_stamp.pl $report > $report_spf");

my $stamp = $report_spf; 
$stamp =~ s/\.spf$/_rel-ab.spf/;

open( 'SPF' , '>' , $stamp ) or die "cant create new spf file $stamp\n";

my %counts = (); ### hash with samples as keys and total read matches as values
my @samples = ();

### read in file once to get sums:
my $lc = 0; # line count
open( 'REPORT' , '<' , $report_spf) or die "cant open kraken-mpa-report output in spf format $report_spf\n";
while( <REPORT> )	{
	
	my @s = split( '	' , $_ );
	
	if ( $lc == 0 )	{
		my @header = (splice( @s , 0 , 7 ) );

		my $i = 0;
		foreach my $s ( @s )	{
				my ($name,$path,$suffix) = fileparse( $s , ("_out.txt") );
				push( @header , $name );
				$counts{$name} = 0;
				push( @samples , $name );
				++$i;
		}
		my $header = join("	" , @header );
		print SPF "$header";
		++$lc;
		next;
	} 

	splice( @s , 0 , 7 );
	
	my $j = 0;
	foreach my $s ( @s )	{
		my $sample = $samples[$j];
		$counts{$sample} += $s;
		++$j;
	}

} close( 'REPORT' );

$lc = 0; # line count
open( 'REPORT' , '<' , $report_spf ) or die "cant open kraken-mpa-report output in spf format $report_spf\n";
while( <REPORT> )	{
	
	my @s = split( '	' , $_ );
	
	if ( $lc == 0 )	{
		++$lc;
		next;
	} 
	
	my @line = splice( @s , 0 , 7 );
	
	my $j = 0;
	foreach my $s ( @s )	{
		my $sample = $samples[$j];
		push( @line , sprintf( "%.5f" , ($s/$counts{$sample})*100) );
		++$j;
	}

	my $line = join( "	" , @line );
	print SPF "$line\n";

} close( 'REPORT' );

### create single log file for all samples by looping through individual log files
open( 'MASTERLOG' , '>' , $log ) or die "cant create MASTERLOG file $log\n";
print MASTERLOG "sample	total_reads	classified	unclassified	classified_percent	unclassified_percent\n";

my $i = 0;

foreach my $err ( @log_files )	{

	### read in log file for kraken jobs:
	open( 'LOG' , '<' , $err ) or die "cant open LOG $err\n";
	chomp( my @lines = <LOG> );
	close( 'LOG' );
	
	my ($name,$path,$suffix) = fileparse( $err , (".fq" , ".fastq") );
	
	### scan 2nd last line of file for info on classified sequences:
	$lines[$#lines-1] =~ m/(\d+) sequences classified \((\S+)%\)/;
	my $classified = $1;
	my $classified_percent = $2;
	
	### scan last line of file for info on unclassified sequences:
	$lines[$#lines] =~ m/(\d+) sequences unclassified \((\S+)%\)/;
	my $unclassified = $1;
	my $unclassified_percent = $2;
	
	my $total = $classified + $unclassified;
	
	print MASTERLOG "$name	$total	$classified	$unclassified	$classified_percent	$unclassified_percent\n";

	### remove individual log files for each vsearch job:
	if ( ! $keep )	{	system( "rm $err" );	}
	++$i;
}

close( 'MASTERLOG' );


sub fastqs2sample {
	
	foreach my $f ( @_ )	{
		my @f_split = ();
		if ( $delimiter eq "_" )	{
			@f_split = split( '_' , basename($f) );
		} elsif ( $delimiter eq "." )	{
			@f_split = split( '\.' , basename($f) );
		}
		my $s = $f_split[0]; #sample
		my $base = basename($f);

		if ( ! exists $s2f{$s} )	{
			my @dummy = ();
			$s2f{$s} = \@dummy;
		}
		push( @{$s2f{$s}} , $f );
	}

	### keep track of #s of SE and PE samples, shouldn't be both types of samples
	my $se_marker = 0;
	my $pe_marker = 0;
	
	### check that if there are 2 fastqs per sample that they are R1 and R2. Also, no more than 2 fastqs per sample.
	foreach my $s ( keys %s2f )	{
		my @f = @{$s2f{$s}};
		my $num_fastq = scalar @f;
		if ( $num_fastq == 1 )	{
			++$se_marker;
			### only 1 fastq for this sample, assuming it's SE and moving on.
		} elsif ( $num_fastq == 2 )	{
			++$pe_marker;
			### assuming this sample has 2 fastqs because they are PE, check that either "_R1." or "_R1_" or "forward" is in one filename and "_R2." or "_R2_" or "reverse" is in the other's filename
			if ( ( basename($f[0]) =~ m/_R1\.|_R1_|forward/ ) and ( basename($f[1]) =~ m/_R2\.|_R2_|reverse/ ) ) {
				### the filenames seem to be PE reads and already ordered as forward and reverse, so continue
			} elsif ( ( basename($f[1]) =~ m/_R1\.|_R1_|forward/ ) and ( basename($f[0]) =~ m/_R2\.|_R2_|reverse/ ) ) 	{
				### filenames are PE, but reorder so they are forward and reverse
				my @tmp = ($f[1] , $f[0] );
				@{$s2f{$s}} = \@tmp;
			} else {
				die "the 2 fastqs for sample $s don't have PE read ids in their names, are you sure these are PE reads? They need _R1. and _R2. (or either _R1_ and _R2_ or forward and reverse) in their names. The filenames are: @f\n";
			}
		} elsif ( ($num_fastq < 1 ) or ( $num_fastq > 2 ) )	{
			die "$num_fastq fastqs for sample $s. This script assumed that fastq prefixes (when delimited by \"$delimiter\") are sample IDs, check your filenames\n";
		}
	}

	if ( ( $se_marker > 0 ) and ( $pe_marker > 0 ) )	{	
		die "Stopping job - $se_marker SE sample(s) and $pe_marker PE sample(s) input. You should run SE and PE samples through this script separately\n";	
	} elsif ( ( $se_marker > 0 ) and ( $paired ) )	{
		die "Stopping job - $se_marker SE sample(s), but '--paired' option used. You should run SE and PE samples through this script separately\n";
	} elsif ( ( $pe_marker > 0 ) and ( ! $paired ) ) 	{
		die "Stopping job - $pe_marker PE sample(s), but '--paired' option was not used.\n";
	}
}

=head1 Name

run_kraken.pl - wrapper to run kraken to do taxonomic classification of metagenomic reads along with the post-processing tools kraken-translate and kraken-mpa-report.  

=head1 USAGE

run_kraken.pl [-log <logfile> -thread <#_CPU_to_use> -o <out_dir> -h -v --keep --preload --quick --min-hits <int> --report <outfile>  --paired --delimiter <.|_>] --db <FASTA> <list of FASTA files>

NOTE: the "kraken" and "kraken-mpa-report" binaries need to be in your PATH.

Four output types will be created:

- a folder containing the output of each individual kraken job (in "kraken_output" by default)
- the output of kraken-mpa-report ("kraken_mpa_report.txt" by default)
- this report converted to STAMP format ("kraken_mpa_report.spf" by default)
- the above STAMP file converted into the relative abundance per sample ("kraken_mpa_report_rel-ab.spf" by default, which is what you would read into STAMP). 

kraken website: 
https://ccb.jhu.edu/software/kraken/

kraken paper: 
Derrick E Wood and Steven L Salzberg
Kraken: ultrafast metagenomic sequence classification using exact alignments
Genome Biology 2014 15:R46
DOI: 10.1186/gb-2014-15-3-r46

=head1 OPTIONS

=over 4

=item B<-h, --help>

Displays the entire help documentation.

=item B<-v, --version>

Displays version number and exits.

=item B<-o, --out_dir <file>>

Output directory for kraken output (default is "kraken_output").

=item B<-thread <# of CPUs>>

Numbers of threads to use (1 by default).

=item B<-keep>

Flag to indicate that temporary log files should not be removed (useful for troubleshooting).

=item B<-log <file>>

the name of the log file - kraken_log.txt by default.
 
=item B<-db, --database <file>> 

Database of genomes to use as a reference (FASTA file).

=item B<--preload > 

kraken flag to preload database into memory for faster computation (off by default).

=item B<--quick > 

kraken flag to indicate quick operation (off by default).

=item B<--min-hits <int> > 

kraken parameter for number of hits required for classification.

=item B<--paired > 

kraken parameter to indicate fastqs are PE.

=item B<--check-names > 

kraken parameter to check that each pair of reads have names that agree with each other. Ignored if --paired is not specified.

=item B<--keep > 

Flag to keep temporary logfiles rather than to delete them (useful for troubleshooting).

=item B<--report <FILE> > 

kraken-mpa-format output file ("kraken_mpa_report.txt" by default)

=item B<-d,--delimiter <.|_>>

Either "." or "_" to use as a delimiter for filenames. Sample IDs are taken to be the first field of each filename. Default: "_".

=back

=head1 AUTHOR

Gavin Douglas <gavin.douglas@dal.ca> 

=cut
