#!/usr/bin/perl

use warnings;
use strict;
use Parallel::ForkManager;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $help; 

my $quality = 0;
my $percent = 0;
my $length = 0;
my $parallel;

my $forward = "ACGCGHNRAACCTTACC";
my $reverse = "TTGYACWCACYGCCCGT";
### (note that the above is the reverse complement of the below primer)
### my $reverse = "ACGGGCRGTGWGTRCAA";

my $out_dir = "filtered_reads";
my $log = "readFilter_log.txt";

my $bbmap_dir = "/usr/local/prg/bbmap";

my $res = GetOptions("out_dir|o=s" => \$out_dir,
		     "log=s"=>\$log,
		     "thread:i"=>\$parallel,
		     "min_quality|q=i" => \$quality,
		     "percent|p=i" => \$percent,
		     "min_length|l=i" => \$length,
		     "help|h"=>\$help,
		     "forward|f=s" => \$forward,
		     "reverse|r=s" => \$reverse,
		     "bbmap|b=s" => \$bbmap_dir,
	  )	or pod2usage(2);

pod2usage(-verbose=>2) if $help;


if ( ( $quality == 0 ) or ( $percent == 0 ) or ( $length == 0 ) )	{
	die "min_quality, percent and min_length are required parameters that need non-zero interger values\nfor help type:	 perl readFilter.pl -h\n";
}


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

system("mkdir -p $out_dir"); ### "-p" makes parent directories as needed

my @files=@ARGV;

pod2usage($0.': You must provide a list of fastq files to be filtered.') unless @files;


my @cmds = ();
my @tmpLog = ();


foreach my $path ( @files )	{
	
	my $file = basename( $path );
	
	my $ext;
	# first check whether filename matches has ".fastq" or ".fq" extension
	if ( $file =~ m/\.fastq$/ )	{
		$ext = "fastq";
	} elsif ( $file =~ m/\.fq$/ )	{
		$ext = "fq";
	} else {	die "file $path does not end in \".fastq\" or \".fq\"\n";	}

	my $outfile = $file;
	$outfile =~ s/\.$ext$/_filtered.$ext/;
	my $outfile_tmp1 = $outfile;
	my $outfile_tmp2 = $outfile;
	my $outfile_tmp3 = $outfile;

	$outfile_tmp1 =~ s/\.$ext$/_TMP1.$ext/;
	$outfile_tmp2 =~ s/\.$ext$/_TMP2.$ext/;
	$outfile_tmp3 =~ s/\.$ext$/_TMP3.$ext/;	

	my $log_tmp = $outfile;
	$log_tmp =~ s/\.$ext$/_TMP_LOG.txt/;

	my $output = $out_dir . "/" . $outfile;
	my $output_tmp1 = $out_dir . "/" . $outfile_tmp1;
	my $output_tmp2 = $out_dir . "/" . $outfile_tmp2;
	my $output_tmp3 = $out_dir . "/" . $outfile_tmp3;
	my $log_tmp_out = $out_dir ."/" . $log_tmp;
	
	my $f_l = length $forward;
	my $r_l = length $reverse;

	my $qFilterCmd = "fastq_quality_filter -v -Q33 -q $quality -p $percent -i $path -o $output_tmp1  >>$log_tmp_out";
	my $lFilterCmd = "$bbmap_dir/bbduk.sh -Xmx1g in=$output_tmp1 outu=$output_tmp2 minlength=$length 2>>$log_tmp_out";
	my $forwardPrimerCmd = "$bbmap_dir/bbduk.sh -Xmx1g in=$output_tmp2 outm=$output_tmp3  restrictleft=$f_l k=$f_l literal=$forward mm=f rcomp=f copyundefined 2>>$log_tmp_out";
	my $reversePrimerCmd = "$bbmap_dir/bbduk.sh -Xmx1g in=$output_tmp3 outm=$output  restrictright=$r_l k=$r_l literal=$reverse mm=f rcomp=f copyundefined 2>>$log_tmp_out";

	my @tmp = ( $qFilterCmd , $lFilterCmd , $forwardPrimerCmd , $reversePrimerCmd , "rm $output_tmp1" , "rm $output_tmp2" , "rm $output_tmp3" );

	push( @cmds , \@tmp );

	push( @tmpLog , "$log_tmp_out,$file" );
}

my $pm = new Parallel::ForkManager($cpu_count);
foreach my $cmds ( @cmds )	{

	$pm->start and next;

	my @c = @{$cmds};

	foreach my $c ( @c )	{
		print STDERR "running: $c\n\n";
		die if system( $c );
	}

	$pm->finish;
}
$pm->wait_all_children;

### parsing logfiles is not paralleled since writing to same file from mutliple jobs can screw up formatting
open( 'LOG' , '>' , $log ) or die "cant create LOG $log\n";
print LOG "file	initial	qFiltered	lFiltered	forwardFiltered	reverseFiltered	final	qFilteredPercent	lFilteredPercent	forwardFilteredPercent	reverseFilteredPercent	finalPercent\n";
foreach my $tmp ( @tmpLog )	{
	&add2log( $tmp );
}
close( 'LOG' );

sub add2log	{

	### take input raw, tmp log file and add it to the cleaned up global log file for all input files
	
	my @s = split( ',' , $_[0] );

	my $tmp = $s[0];
	my $name = $s[1];

	my @inCount = ();
	my @outCount = ();

	my $resultCount = 0;

	open( 'TMP' , '<' , $tmp ) or die "cant open TMP logfile $tmp\n";
	while( <TMP> )	{
		
		my @split = split( '\s+' , $_ );
		
		if ( ! exists $split[0] )	{	next	}

		if ( $split[0] eq "Input:" )	{
			push( @inCount , $split[1] );
		} elsif ( ( $split[0] eq "Contaminants:" ) or ( $split[0] eq "Output:" ) or ( $split[0] eq "Result:" ) )	{
			if ( $split[0] eq "Result:" )	{
				if ( $resultCount > 0 )	{
					next;
				} else {
					++$resultCount;
				}
			} else {}
			push( @outCount , $split[1] );
		} else {}
		
	} close( 'TMP' );

	my $initial = $inCount[0];
	my $qFiltered = $inCount[0] - $inCount[1];
	my $lFiltered = $inCount[1] - $inCount[2];
	my $forwardFiltered = $inCount[2] - $inCount[3];
	my $reverseFiltered = $inCount[3] - $outCount[3];
	my $final = $outCount[3];

	### note that percents are all based on initial count!
	my $qPercent = sprintf( "%.1f" , ($qFiltered / $initial)*100  );
	my $lPercent = sprintf( "%.1f" , ($lFiltered / $initial)*100  );
	my $forwardPercent = sprintf( "%.1f" , ($forwardFiltered / $initial)*100  );
	my $reversePercent = sprintf( "%.1f" , ($reverseFiltered / $initial)*100  );
	my $finalPercent = sprintf( "%.1f" , ($final / $initial)*100 );

	print LOG "$name	$initial	$qFiltered	$lFiltered	$forwardFiltered	$reverseFiltered	$final	$qPercent	$lPercent	$forwardPercent	$reversePercent	$finalPercent\n";

	system( "rm $tmp" );
}

__END__

=head1 Name

readFilter.pl - wrapper to filter reads by quality with fastx and then by total length with bbmap. Reads without matches to the forward and reverse primers are then removed with bbmap. 

=head1 USAGE

readFilter.pl [-f <oligo> -r <oligo> -bbmap <directory> -log <logfile> -thread <#_CPU_to_use> -o <out_dir> -h] -q <min_quality> -p <min_percent_sites_with_q> -l <min_length>  <list of fastq files>


Examples:

# remove all reads that do not have a quality score of 30 at least 90% of bases. Then remove all reads that are less than 400bp long.

readFilter.pl -q 30 -p 90 -l 400 *.fastq


# thread with 2 CPUs, remove all reads that do not have a quality score of 30 at least 90% of bases. Then remove all reads that are less than 400bp long.

readFilter.pl -thread 2 -q 30 -p 90 -l 400 *.fastq


# thread on all available CPUs, output into "filtered_reads", write log output to "filtered.log", min quality of 20, min percentage of bases with that quality of 90%, min length of 350 bases.

readFilter.pl -thread -o filtered_reads -log filtered.log -q 20 -p 80 -l 350 *.fastq


=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

Output directory for filtered fastq files. Default is "filtered_reads".

=item B<-thread [<# of CPUs>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-log <file>>

The location to write the log file.
 
=item B<-h, --help>

Displays the entire help documentation.

=item B<-q, --min_quality>

Minimum base quality.

=item B<-p, --percent>

Minimum percent of bases per read that pass quality cut-off

=item B<-l, --min_length>

Minimum read length.

=item B<-f, --forward>

Forward primer to match at beginning of all reads (IUPAC format, default: ACGCGHNRAACCTTACC).

=item B<-r, --reverse>

Reverse primer to match at end of all reads (IUPAC format, default: TTGYACWCACYGCCCGT, which is the reverse complement of the primer ACGGGCRGTGWGTRCAA).

=item B<-b, --bbmap>

bbmap directory containing sh files (default: /usr/local/prg/bbmap). 

=back

=head1 DESCRIPTION

B<readFilter.pl> This script automatically filters multiple fastqs by quality and length.

The script allows the use of multiple threads. 

By default, log output is written to "readFilter_log.txt".

bbmap is hard coded into this script, so this will have to changed on a different system (see "--bbmap" option). Also, FASTX-Toolkit needs to be installed and in the user's $PATH.


# software websites:
http://sourceforge.net/projects/bbmap/
http://hannonlab.cshl.edu/fastx_toolkit/


=head1 AUTHOR

Gavin Douglas <gavin.douglas@dal.ca> (based on structure by Morgan Langille)

=cut

