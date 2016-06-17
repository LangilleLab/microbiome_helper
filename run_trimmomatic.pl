.#!/usr/bin/perl

use warnings;
use strict;

use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $help; 
my $version_marker;
my $version = "1.0";

### default cut-offs
my $leading_min_qual = 5;
my $trailing_min_qual = 5;
my $window_size = 4;
my $required_qual = 15;
my $min_length = 70;
 
### by default will run on 1 CPU
my $parallel = 1;

### will not overwrite by default
my $force;

### will delete individual log files by default
my $keep;

### default paths:
my $out_dir = "trimmomatic_filtered";
my $log = "trimmomatic_tabular_log.txt";
my $trimmomatic_jar = "/usr/local/prg/Trimmomatic-0.36/trimmomatic-0.36.jar";

### delimiter used to identify sample name (either "_" or ".")
my $delimiter = "_";

my $res = GetOptions(
		       	"out_dir|o=s" => \$out_dir,
		         "log=s"=>\$log,
		         "thread:i"=>\$parallel,
		         "required_quality|r=i" => \$required_qual,
	         	 "leading_quality|l=i" => \$leading_min_qual,
			 "trailing_quality|t=i" => \$trailing_min_qual,
			 "window_size|w=i" => \$window_size,
			 "min_length|m=i" => \$min_length,
			 "jar|j=s" => \$trimmomatic_jar,
			 "force|f" => \$force,
			 "keep|k" => \$keep,
			 "delimiter|d=s" =>\$delimiter,
		         "help|h"=>\$help,
		         "version|v" => \$version_marker,
	  )	or pod2usage(2);

pod2usage(-verbose=>2) if $help;

if ( $version_marker )	{	print "version $version\n";	exit	}

if ( ! -e $trimmomatic_jar )	{	die "$trimmomatic_jar does not exist, you may need to set the --jar flag\n";	}

if ( ( $delimiter ne "_" ) and ( $delimiter ne "." ) )	{	die "delimiter $delimiter needs to be \"_\" or \".\"\n";	}

my @files = @ARGV; 

pod2usage($0.': You must provide a list of fastq files to be trimmed. Use the -h option for more info.') unless @files;

### first check that all input fastqs exist:
foreach my $f ( @files )	{ 
	if ( ! -e $f ) { die "Stopping: file $f not found\n";	}
}

### figure out fastqs for each sample (important because there could be paired-end fastqs input).

my %s2f = (); ###sample2fastqs
&fastqs2sample( @files ); ### this subroutine will add fastqs to above hash with arrays as values and samples as keys

system("mkdir -p $out_dir"); ### "-p" makes parent directories as needed

my @cmds = ();
my @tmpLog = ();

foreach my $s ( keys %s2f )	{	### loop over sample names
	
	print STDERR "\nstarting on sample $s\n";
	
	my @f = @{$s2f{$s}}; ### files for sample $s
	my $file1 = $f[0];

	# check whether filename matches has ".fastq" or ".fq" extension
	my $ext;
	if ( $file1 =~ m/\.fastq$/ )	{
		$ext = "fastq";
	} elsif ( $file1 =~ m/\.fq$/ )	{
		$ext = "fq";
	} else {	die "file $file1 does not end in \".fastq\" or \".fq\"\n";	}

	my $steps = "LEADING:$leading_min_qual TRAILING:$trailing_min_qual SLIDINGWINDOW:$window_size:$required_qual MINLEN:$min_length 2>$out_dir/$s"."_tmp_trim.log";
	my $opt = "-threads $parallel -phred33";
	my $jar_cmd = "java -jar $trimmomatic_jar";

	push( @tmpLog , "$out_dir/$s"."_tmp_trim.log,$s" );

	if ( scalar @f == 1 )	{	### SE reads
		
		my $outfile = $out_dir . "/" . basename($file1);
		$outfile =~ s/\.$ext$/_trimmed.$ext/;
		
		if ( (-e $outfile ) and ( ! defined $force ) )	{	die "stopping: output file $outfile already exists, use the --force flag if you want to overwrite previous files\n";	}
		
		# run command:
		print STDERR "running: $jar_cmd SE $opt $file1 $outfile $steps\n";
		system( "$jar_cmd SE $opt $file1 $outfile $steps");
		
	} elsif ( scalar @f == 2 )	{
		
		my $file2 = $f[1];
		
		my $outfile_paired1 = $out_dir . "/" . basename($file1);
		$outfile_paired1  =~ s/\.$ext$/_trimmed_paired.$ext/;
		
		my $outfile_paired2 = $out_dir . "/" . basename($file2);
		$outfile_paired2  =~ s/\.$ext$/_trimmed_paired.$ext/;
		
		my $outfile_unpaired1 = $out_dir . "/" . basename($file1);
		$outfile_unpaired1  =~ s/\.$ext$/_trimmed_unpaired.$ext/;
		
		my $outfile_unpaired2 = $out_dir . "/" . basename($file2);
		$outfile_unpaired2  =~ s/\.$ext$/_trimmed_unpaired.$ext/;
		
		if ( ( (-e $outfile_paired1 ) or ( -e $outfile_unpaired1 ) or ( -e $outfile_paired2 ) or ( -e $outfile_unpaired2 ) ) and ( ! defined $force ) )	{	die "stopping: output files for sample $s already exist, use the --force flag if you want to overwrite previous files\n";	}
		
		#run command:
		print STDERR "running: $jar_cmd PE $opt $file1 $file2 $outfile_paired1 $outfile_unpaired1 $outfile_paired2 $outfile_unpaired2 $steps\n";
		system( "$jar_cmd PE $opt $file1 $file2 $outfile_paired1 $outfile_unpaired1 $outfile_paired2 $outfile_unpaired2 $steps");
		
	}	
		
}		
		
### parsing logfiles is not paralleled since writing to same file from mutliple jobs can screw up formatting
open( 'LOG' , '>' , $log ) or die "cant create LOG $log\n";

my $h_marker = 0;

foreach my $tmp ( @tmpLog )	{
	&add2log( $tmp );
}
close( 'LOG' );

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

	if ( ( $se_marker > 0 ) and ( $pe_marker > 0 ) )	{	die "$se_marker SE sample(s) and $pe_marker PE sample(s) input. You should run SE and PE samples through this script separately\n";	}
	
}

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

		if ( $split[0] eq "Input" )	{

			if ( $_ =~ m/Input Read Pairs: (\d+) Both Surviving: (\d+) \((\S+)%\) Forward Only Surviving: (\d+) \((\S+)%\) Reverse Only Surviving: (\d+) \((\S+)%\) Dropped: (\d+) \((\S+)%\)/ )	{
	
				if ( $h_marker == 0 )	{	
					++$h_marker;
					print LOG "sample	input	both_passed	forward_only_passed	reverse_only_passed	dropped	both_passed_percent	forward_only_passed_percent	reverse_only_passed_percent	dropped_percent\n";
				} else {}
	
				print LOG "$name	$1	$2	$4	$6	$8	$3	$5	$7	$9\n";

			} elsif ( $_ =~ m/Input Reads: (\d+) Surviving: (\d+) \((\S+)%\) Dropped: (\d+) \((\S+)%\)/)	{
			
				if ( $h_marker == 0 )	{
					++$h_marker;
					print LOG "sample	input	passed	passed_percent\n";
				}

				print LOG "$name	$1		$2	$3\n";

			} else {}

		}

	} close ('TMP');
	if ( ! $keep )	{	system( "rm $tmp" )	}
}

__END__

=head1 Name

run_trimmomatic.pl - wrapper to trim reads with Trimmomatic.

=head1 USAGE

run_trimmomatic.pl  [-l <int> -t <int> -r <int> -w <int> -m <int> --jar <Path to Trimmomatic jarfile> --log <logfile> --thread <#_CPU_to_use> -o <out_dir> -h -v] <list of fastq files>

Note: all fastq files containing "forward", "reverse", "_R1.", "_R1_", "_R2." or "_R2_" will be interpreted as paired-end reads, otherwise they will be assumed to be single ended. The sample ID is also assumed to be the prefix of each filename when delimited by "_".

Examples:

# Run with all defaults:

run_trimmomatic.pl *.fastq

# Trim leading and trailing bases that have qualities less than 5 (-l and -t options), cut at windows of size 4 (-w option) that have an average quality less than 15 (-r option) and throw out reads that are less than 70 nucleotides after trimming (-m option). Run trimmomatic with the specified jarfile (-j option) using 20 CPUs (--thread option) and overwrite previous output files (-f option):

run_trimmomatic.pl -l 5 -t 5 -r 15 -w 4 -m 70 -j /usr/local/prg/Trimmomatic-0.36/trimmomatic-0.36.jar -f --thread 20 *.fastq


=head1 OPTIONS

=over 4

=item B<-h, --help>

Displays the entire help documentation.

=item B<-v, --version>

Displays script version and exits.

=item B<-o, --out_dir <file>>

Output directory for filtered fastq files (default: "trimmomatic_filtered").

=item B<--thread <# of CPUs>>

Number of CPUs to thread each job on (default: 1).

=item B<--log <file>>

Name of log file (default: trimmomatic_tabular_log.txt).
 
=item B<-l,--leading_quality <int>>

The minimum quality for leading bases to be kept as required by Trimmomatic's LEADING command (default: 5).

=item B<-t,--trailing_quality <int>>

The minimum quality for trailing bases to be kept as required by Trimmomatic's TRAILING command (default: 5).
 
=item B<-r,--required_quality <int>>

Average quality required by Trimmomatic's SLIDINGWINDOW command (default: 15).

=item B<-w,--window_size <int>>

Window sizes of sliding windows required by Trimmomatic's SLIDINGWINDOW command (default: 4).

=item B<-m,--min_length <int>>

Min length of reads to be retained as required by Trimmomatic's MINLEN command (default: 70).

=item B<-j,--jar <jarfile>>

Path to Trimmomatic jarfile (default: /usr/local/prg/Trimmomatic-0.36/trimmomatic-0.36.jar).

=item B<-f,--force >

Flag to indicate previous files should be overwritten.

=item B<-k, --keep >
	
Flag to indicate that individual log files should not be deleted (helpful for troubleshooting).
	
=item B<-d,--delimiter <.|_>>

Either "." or "_" to use as a delimiter for filenames. Sample IDs are taken to be the first field of each filename. Default: "_".

=back

=head1 DESCRIPTION

B<run_trimmomatic.pl> This script wraps Trimmomatic to trim reads.

The script allows the use of multiple threads. 

By default, log output is written to "trimmomatic_tabular_log.txt".

As stated above, note that all fastq files containing "forward" , "_R1.", "_R1_", "reverse", "_R2." or "_R2_" will be interpreted as paired-end reads, otherwise they will be assumed to be single ended. 

Also, the sample ID is assumed to be the first field of each filename when delimited by "_".


Trimmomatic citation:
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Trimmomatic website: 
http://www.usadellab.org/cms/index.php?page=trimmomatic


=head1 AUTHOR

Gavin Douglas <gavin.douglas@dal.ca> 

=cut

