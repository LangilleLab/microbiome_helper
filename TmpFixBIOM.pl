#!/usr/bin/perl

use warnings;
use strict;

if ( scalar @ARGV != 1 )	{	die "only input argument needs to be BIOM file in TSV format\n";	}

### fix taxonomies that include "uncultured", "Incertae Sedis", "environmental sample"

my $lc = 0;
open ( 'BIOM' , '<', $ARGV[0] ) or die "cant open BIOM $ARGV[0]\n";
while( <BIOM> )	{
	my @s = split( '[\t\n]' , $_ );
	my $first = substr( $_ , 0 , 1 );
	if ( $first eq "#" )	{	print "$_";	next	}
	
	my $otu = $s[0];
	my @l = split( ';\s' , $s[$#s] );
	my $last;
	for ( my $i = 0; $i < scalar @l ; ++$i )	{
		if ( $l[$i] =~ m/\s/g )	{	$l[$i] =~ s/ /_/g;	}
		if ( ( $l[$i] =~ m/uncultured/g ) or ( $l[$i] =~ m/environmental_sample/g ) )	{	$l[$i] = $l[$i] . "_".$otu }
		if ( $l[$i] =~ m/Incertae_Sedis/g )	{	$l[$i] = $l[$i] . "_".$last } 
		$last = (split( '__' , $l[$i]))[1];
	}
	$s[$#s] = join( "; " , @l );
	my $line = join( "	" , @s );
	print "$line\n"; 
} close( 'BIOM' );
