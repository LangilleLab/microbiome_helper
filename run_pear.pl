#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use List::Util qw(min max sum);


my ($parallel,$help);
my $out_dir='./';
my $full_log='pear_full_log.txt';
my $summary_log='pear_summary_log.txt';
my $gzip_output;

my $res = GetOptions("out_dir=s" => \$out_dir,
		     "parallel:i"=>\$parallel,
		     "full_log=s"=>\$full_log,
		     "summary_log=s"=>\$summary_log,
		     "gzip_output"=>\$gzip_output,
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
    #attempt different naming scheme
    }elsif($file_name =~ /(.+)_([1|2])/){
	$paired_files{$1}[$2-1]=$file;
    }else{
	warn "Input file does not contain '_R1_' or '_R2_' in name: $file";
    }
}

#clear the output log (and make sure it is writable)
#open(my $FULL_LOG,'>',$full_log) || die "Can't write to log file: $full_log";
#close($FULL_LOG);


foreach my $name (sort keys %paired_files){
    unless(defined($paired_files{$name}[0]) && defined($paired_files{$name}[1])){
	warn "Couldn't find matching paired end files for file starting with: $name";
	next;
    }
    my $out_file=$out_dir.'/'.$name;
    #check if this has already been done
    my $assembled_out_file=$out_file.'.assembled.fastq';
    if (-e $assembled_out_file || -e $assembled_out_file.'.gz'){
	print "Skipping this sample because output file already exists: $assembled_out_file\n";
	next;
    }
    my $cmd="pear -f $paired_files{$name}[0] -r $paired_files{$name}[1] -j $cpu_count -o $out_file >>$full_log";
    print $cmd,"\n";
    die if system($cmd);

    #compress output files (if the flag is set)
    if($gzip_output){
	my $gzip_cmd="pigz -p $cpu_count -f $out_file".'*';
	print $gzip_cmd,"\n";
	die if system($gzip_cmd);
    }
}

print "Creating PEAR summary log at: $summary_log \n";
my $min_assembled=create_summary_log($full_log,$summary_log);

if($min_assembled < 90){
    print "Finished! Warning!! one or more samples were less than 90% assembled! You should manually inspect the log file: $summary_log \n";
}else{
    print "Finished! All samples assembled at 90% or greater. For more details you can check manually inspect the log file: $summary_log \n";
}

sub mean {
    return sum(@_)/@_;
}

sub create_summary_log{

    my $full_log=shift;
    my $summary_log=shift;
    open(my $FULL_LOG,'<',$full_log) || die "Can't read log file: $full_log";

    open(my $SUMMARY_LOG,'>',$summary_log) || die "Can't create summary log file for writing: $summary_log";
    my @samples;
    
    while (<$FULL_LOG>) {
	chomp;
	if (/Assembled reads/) {
	    my $assembled_string=$_;
	    my $discarded_string=<$FULL_LOG>;
	    my $unassembled_string=<$FULL_LOG>;
	    my $assembled_file_string=<$FULL_LOG>;
	    my ($assembled_percent) = $assembled_string =~ /(\d+\.\d+)\%/;
	    my ($discarded_percent) = $discarded_string =~ /(\d+\.\d+)\%/;
	    my ($unassembled_percent) = $unassembled_string =~ /(\d+\.\d+)\%/;
	    my ($assembled_file) = $assembled_file_string =~ /([\w|\-|\.]+)\.assembled/;
	    push (@samples, [$assembled_file,$assembled_percent,$discarded_percent,$unassembled_percent]);
	}
    }

    #Add min, mean, and max as first three lines of output
    unshift @samples,['Max',sprintf("%.3f",max(map{$_->[1]}@samples)),sprintf("%.3f",max(map{$_->[2]}@samples)),sprintf("%.3f",max(map{$_->[3]}@samples))];
    unshift @samples,['Mean',sprintf("%.3f",mean(map{$_->[1]}@samples)),sprintf("%.3f",mean(map{$_->[2]}@samples)),sprintf("%.3f",mean(map{$_->[3]}@samples))];
    unshift @samples,['Min',sprintf("%.3f",min(map{$_->[1]}@samples)),sprintf("%.3f",min(map{$_->[2]}@samples)),sprintf("%.3f",min(map{$_->[3]}@samples))];

    #print header
    print $SUMMARY_LOG join("\t","ID","Assembled","Discarded","Unassembled"),"\n";

    #print out all the data
    foreach my $sample (@samples) {
	print $SUMMARY_LOG join("\t",@$sample),"\n";
    }
    return sprintf("%.3f",min(map{$_->[1]}@samples))
}


__END__

=head1 Name

run_pear.pl - A simple wrapper for PEAR to stich paired-end reads

=head1 USAGE

run_pear.pl [-p [<# proc>] -o <out_dir> -h] <list of fastq files>

E.g.

#Note: Files must have "_R1_" and "_R2_" within the file name (or secondarily "_1" and "_2")

run_pear.pl sample1_R1_001.fastq sample1_R2_001.fastq sample2_R1_001.fastq sample2_R2_001.fastq

#Shorter way to do the same thing

run_pear.pl *.fastq

#Specify alternate location for output files (instead of default current directory)

run_pear.pl -o stitched_reads *.fastq

#Run in parallel and use all CPUs

run_pear.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

run_pear.pl *.fastq -p 2

#Turn off gzip compression of output files

run_pear.pl -g *.fastq

=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all PEAR output files.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-g, --gzip_output>

Gzip the PEAR output files.

=item B<-f, --full_log <file>>

The location to write the PEAR full log file. Default is "pear_full_log.txt"

=item B<-s, --summary_log <file>>

The location to write teh PEAR summary log file. Default is "pear_summary_log.txt"

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<run_pear.pl> This script allows for more automated running of the PEAR program on multiple fastq files. PEAR is used to stitch (or assemble) paired end reads together. The assumption is made that the paired end files have the same name with the forward reads being indicated by "_R1_" and the reverse being "_R2_". If file names are not found matching these then an simpler label is attempted ("_1" and "_2"). 

The script allows the use of multiple threads. 

This script also captures the output statistics from PEAR and outputs them to a single "pear_full_log.txt"(by default). It also parses this and simplifies the output into "pear_summary_log.txt" (by default). 
 
By default, output files from PEAR are gzipped to save on space. 

Before use make sure you have installed the "pear" program so it is accesible from your PATH.

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

