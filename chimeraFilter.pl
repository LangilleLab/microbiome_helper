#!/usr/bin/perl

use warnings;
use strict;
use Parallel::ForkManager;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $help; 

my $parallel;
my $db ;
my $out_dir = "non_chimeras";
my $log = "chimeraFilter_log.txt";
my $type;
my $mindiv = 1.5;
my $minh = 0.2;

my $res = GetOptions("type:i"=> \$type,
		     "out_dir|o=s" => \$out_dir,
		     "database|db=s"=>\$db,
		     "log=s"=>\$log,
		     "thread:i"=>\$parallel,
		     "mindiv=f"=>\$mindiv,
		     "minh=f"=>\$minh,
		     "help|h"=>\$help,
	  )	or pod2usage(2);

pod2usage(-verbose=>2) if $help;

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

if ( ( ! defined $type ) or (  ( $type != 0 ) and ( $type != 1 ) ) )	{	die "output type \"-type\" needs to be either 1 (only clear non-chimeras) or 0 (any non-chimera, which includes those which are borderline)\n";	}

if ( -e $db )	{	} else { die "database file $db not found\n";	}

my @files=@ARGV;

pod2usage($0.': You must provide a list of fasta files to be filtered.') unless @files;
if ( ! defined $db )	{	die "need database file\n";	}

open( 'MASTERLOG' , '>' , $log ) or die "cant create MASTERLOG file $log\n";
print MASTERLOG "file	chimeraCalls	unclearCalls	nonChimeraCalls	chimeraCallsPercent	unclearCallsPercent	nonChimeraCallsPercent\n";

system( "mkdir -p $out_dir");

my $pm = new Parallel::ForkManager($cpu_count);

foreach my $f ( @files )	{
		
	$pm->start and next;
	
	my $base = basename( $f );
	my $err = $out_dir ."/".$base . ".LOG";
	
	my @baseSplit = split( '\.' , $base );

	my $ext = pop @baseSplit;
	my $out = join( "." , @baseSplit ) . ".nonchimera.".$ext;

	$out = "$out_dir" ."/".$out;
	
	my $uchimeOutTmp = $out . ".uchimeTMP.txt";
			
	print STDERR "usearch61 -uchime_ref $f -db $db  -minh $minh -mindiv $mindiv -threads 1 -strand plus -uchimeout $uchimeOutTmp 2>$err\n";
	
	system( "usearch61 -uchime_ref $f -db $db  -minh $minh -mindiv $mindiv -threads 1 -strand plus -uchimeout $uchimeOutTmp 2>$err >>uchime_tmp.txt" );
		
	my %calls = ();
	my %counts = ();

	$counts{"Y"} = 0;
	$counts{"N"} = 0;
	$counts{"?"} = 0;

	open( 'UCHIME' , '<' , $uchimeOutTmp ) or die "cant open UCHIME $uchimeOutTmp\n";
	while (<UCHIME> )	{
		my @split = split( '\s+' , $_ );
		$calls{$split[1]} = $split[$#split];
		$counts{$split[$#split]} += 1;
	} close( 'UCHIME' );


	open( 'OUT' , '>' , $out ) or die "cant create OUT file $out\n";
	open ( 'IN' , '<' , $f ) or die "cant open IN file $f\n";

	my $next = 0;

	while ( <IN> )	{
		
		my @split = split( '\s+' , $_ );
		if ( ! exists $split[0] )	{	next	};
	
		my $first = substr( $_ , 0 , 1 );

		if ( $first eq ">" )	{
		
			my $name = substr( $split[0] , 1 );
			
			if ( ! exists $calls{$name} )	{	die "sequence $name not found in uchime output. Did uchime run properly?\n";	}		
	
			if ( $calls{$name} eq "N" )	{
				print OUT "$_";	
				$next = 1;	
			} elsif ( ( $calls{$name} eq "?" ) and ( $type == 0 ) )	{
				print OUT "$_";
				$next = 1;	
			} else {
				$next = 0;
			}
		} elsif ( $next == 1 )	{
			print OUT "$_";
		} else {}
	} close ('IN' );
	close( 'OUT' );

	my $chimeras = $counts{"Y"};
	my $non = $counts{"N"};
	my $unclear = $counts{"?"};

	my $total = $chimeras + $non + $unclear;

	if ( $total == 0 )	{	die "no sequences in uchime output file... stopping\n";	}

	my $perChimera = ( $chimeras / $total) * 100;
	my $perNon = ( $non / $total ) * 100;
	my $perUnclear = ( $unclear / $total) * 100;
	
	$perChimera = sprintf ( "%.1f" , $perChimera );
	$perNon = sprintf ( "%.1f" , $perNon );
	$perUnclear = sprintf ( "%.1f" , $perUnclear );

	print MASTERLOG "$base	$chimeras	$unclear	$non	$perChimera	$perUnclear	$perNon\n";

	system( "rm $uchimeOutTmp" );
	system( "rm $err" );
	
	$pm->finish;
}

$pm->wait_all_children;

system( "rm uchime_tmp.txt");

=head1 Name

chimeraFilter.pl - wrapper to filter out chimeric reads from fasta files (using usearch 6.1, specifically uchime).


=head1 USAGE

chimeraFilter.pl [-log <logfile> -thread <#_CPU_to_use> -o <out_dir> -minh <minimum chimera score> -mindiv <minimum divergence> -h]  -type <0 or 1> -db <database of known 16S genes> <list of fasta files>


NOTE: currently the binary "usearch61" (USEARCH 6.1) needs to be in your path. This will be fixed soon!

=head1 OPTIONS

=over 4

=item B<-type <0 or 1>>

Non-chimeric output type, either only sequences that are clearly non-chimeric (1),

or

all sequences that are not called as chimeric ( 0 - includes borderline sequences, "?" in uchime output).


=item B<-mindiv <number>>

Min % divergence between query and target sequence (default 1.5, note that this differs from the uchime default of 0.8).

=item B<-minh <number>>

Min score to be called as chimeric (default 0.2, note that this differs from the uchime default of 0.28).

=item B<-o, --out_dir <file>>

Output directory for filtered fastq files. Default is "non_chimeras".

=item B<-thread [<# of CPUs>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-log <file>>

The location to write the log file.
 
=item B<-h, --help>

Displays the entire help documentation.

=item B<-db, --database> 

Database of 16S sequences to use as a reference (.udb or fasta file).

=back

=head1 AUTHOR

Gavin Douglas <gavin.douglas@dal.ca> (based on structure by Morgan Langille)

=cut
