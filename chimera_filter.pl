#!/usr/bin/perl

use warnings;
use strict;

use Parallel::ForkManager;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $help; 
my $version = "2.0";
my $version_marker;

my $parallel;
my $db ;
my $out_dir = "non_chimeras";
my $log = "chimera_filter_log.txt";
my $type;
my $mindiv = 1.5;
my $minh = 0.2; 
my $keep;  
my @log_files = (); ### array to contain all log files, which will be looped over at the end

my $res = GetOptions("type:i"=> \$type,
		     "out_dir|o=s" => \$out_dir,
		     "database|db=s"=>\$db,
		     "log=s"=>\$log,
		     "thread:i"=>\$parallel,
		     "mindiv=f"=>\$mindiv,
		     "minh=f"=>\$minh,
		     "keep"=>\$keep,
		     "help|h"=>\$help,
		     "version|v"=>\$version_marker,
	  )	or pod2usage(2);

pod2usage(-verbose=>2) if $help;

if ( $version_marker )	{	print "version $version\n";	exit	}

my $cpu_count=1;
#if the option is set
if(defined($parallel)){
    #option is set but with no value then use the max number of proccessors
    if($parallel ==0){
		#load this module dynamically
		eval("use Sys::CPU;");
		$cpu_count=Sys::CPU::cpu_count();
    } else {
		$cpu_count=$parallel; 
    }
}

### the "2>&1" syntax means output STDERR to same output as STDOUT (i.e. same output as file descripter 1).
if ( index( `vsearch 2>&1` , "torognes/vsearch" ) == -1 )      { die "Stopping job, because \"vsearch\" is not in your path.\n";	} 	

if ( ( ! defined $type ) or (  ( $type != 0 ) and ( $type != 1 ) ) )	{	die "output type \"-type\" needs to be either 1 (only clear non-chimeras) or 0 (any non-chimera, which includes those which are borderline)\n\nType \"chimera_filter.pl -h\" for help.\n\n";	}

if ( ! defined $db )	{	die "need database file\n";	}
if ( ! -e $db )	{ die "specified database $db does not exist\n";	}

my @files = @ARGV;
pod2usage($0.': You must provide a list of fasta files to be filtered.') unless @files;

system( "mkdir -p $out_dir");

my $pm = new Parallel::ForkManager($cpu_count);

foreach my $f ( @files )	{
	
	my $base = basename( $f );
	
	### store log files to loop through at the end:
	my $err = $out_dir ."/".$base . ".LOG";
	push( @log_files , $err );
	
	$pm->start and next;
	
	my @baseSplit = split( '\.' , $base );
	
	my $ext = pop @baseSplit;
	my $out = join( "." , @baseSplit ) . ".nonchimera.".$ext;
	my $unclear_out = join( "." , @baseSplit ) . ".unclear.".$ext;
	
	$out = "$out_dir" ."/".$out;
	$unclear_out = "$out_dir" ."/". $unclear_out;
	
	if ( $type == 0 )	{	
		
		print STDERR "vsearch --uchime_ref $f --db $db --minh $minh --mindiv $mindiv --threads 1 --nonchimeras $out --borderline $unclear_out >> tmp_vsearch_out.txt  2> $err\n"; 
		system( "vsearch --uchime_ref $f --db $db --minh $minh --mindiv $mindiv --threads 1 --nonchimeras $out --borderline $unclear_out >> tmp_vsearch_out.txt  2> $err" );
		
		my $any_non_chimera = "$out_dir" . "/" . join( "." , @baseSplit ) . ".nonchimera_and_unclear." . $ext;		
		
		### combine clear non-chimeras and borderline sequences:
		print STDERR "cat $out $unclear_out > $any_non_chimera\n";
		system( "cat $out $unclear_out > $any_non_chimera" );	
		
		if ( ! $keep )	{
			### remove temp out files:
			print STDERR "rm $out\n";
			system( "rm $out" );
			print STDERR "rm $unclear_out\n";
			system( "rm $unclear_out" );
		}
		
	} elsif ( $type == 1 )	{
	
		print STDERR "vsearch --uchime_ref $f --db $db --minh $minh --mindiv $mindiv --threads 1 --nonchimeras $out >> tmp_vsearch_out.txt  2> $err\n"; 
		system( "vsearch --uchime_ref $f --db $db --minh $minh --mindiv $mindiv --threads 1 --nonchimeras $out >> tmp_vsearch_out.txt  2> $err" );

	}
	
	$pm->finish;
}

$pm->wait_all_children;

### create single log file for all samples by looping through individual log files
open( 'MASTERLOG' , '>' , $log ) or die "cant create MASTERLOG file $log\n";
print MASTERLOG "file	chimeraCalls	unclearCalls	nonChimeraCalls	chimeraCallsPercent	unclearCallsPercent	nonChimeraCallsPercent\n";

my $i = 0;

foreach my $err ( @log_files )	{

	### read in log file for vsearch job:
	open( 'LOG' , '<' , $err ) or die "cant open LOG $err\n";
	chomp( my @lines = <LOG> );
	close( 'LOG' );
	
	my $base = basename( $files[$i] );
	
	### scan 2nd last line of file for info on chimeras and non-chimeras:
	$lines[$#lines-1] =~ m/(\d+) \((\S+)%\) chimeras, (\d+) \((\S+)%\) non-chimeras,/;
	my $chimeras = $1;
	my $perChimera = $2;
	my $non = $3;
	my $perNon = $4;
	
	### scan last line of file for info on borderline (i.e. unclear) sequences:
	$lines[$#lines] =~ m/and (\d+) \((\S+)%\) borderline sequences in (\d+) total sequences\./;
	my $unclear = $1;
	my $perUnclear = $2;
	
	print MASTERLOG "$base	$chimeras	$unclear	$non	$perChimera	$perUnclear	$perNon\n";

	### remove individual log files for each vsearch job:
	if ( ! $keep )	{	system( "rm $err" );	}
	++$i;
}

close( 'MASTERLOG' );
system( "rm tmp_vsearch_out.txt" );

=head1 Name

chimera_filter.pl - wrapper to filter out chimeric reads from fasta files (with vsearch, which uses the uchime algorithm).


=head1 USAGE

chimera_filter.pl [-log <logfile> -thread <#_CPU_to_use> -o <out_dir> -minh <minimum chimera score> -mindiv <minimum divergence> -h -v --keep]  -type <0 or 1> -db <FASTA> <list of FASTA files>


NOTE: the "vsearch" binary needs to be in your PATH.

vsearch GitHub page: https://github.com/torognes/vsearch/  

uchime paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3150044/

uchime algorithm: http://drive5.com/usearch/manual/uchime_algo.html


=head1 OPTIONS

=over 4

=item B<-h, --help>

Displays the entire help documentation.

=item B<-v, --version>

Displays version number and exits.

=item B<-type <0 or 1>>

Non-chimeric output type, either only sequences that are clearly non-chimeric (1),

or

all sequences that are not called as chimeric ( 0 - includes borderline sequences).


=item B<-mindiv <float>>

Min % divergence between query and target sequence (default 1.5, note that this differs from the vsearch default of 0.8).

=item B<-minh <float>>

Min score to be called as chimeric (default 0.2, note that this differs from the vsearch default of 0.28).

=item B<-o, --out_dir <file>>

Output directory for filtered fastq files. Default is "non_chimeras".

=item B<-thread <# of CPUs>>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-keep>

Flag to indicate that temporary log files should not be removed (useful for troubleshooting). Also, will prevent the "nonchimera" and "unclear" specific fastas from being removed when type == 0. 

=item B<-log <file>>

The location to write the log file.
 
=item B<-db, --database <file>> 

Database of 16S sequences to use as a reference (FASTA file).

=back

=head1 AUTHOR

Gavin Douglas <gavin.douglas@dal.ca> (based on structure by Morgan Langille)

=cut
