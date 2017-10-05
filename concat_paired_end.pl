#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;


my ($gzip,$write_only,$parallel,$no_R_match,$full_name,$help);
my $out_dir='./';
my $res = GetOptions("out_dir=s" => \$out_dir,
		     "parallel:i"=>\$parallel,
		     "gzip"=>\$gzip,
		     "write_only"=>\$write_only,
                     "no_R_match"=>\$no_R_match,
		     "full_name"=>\$full_name,
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
    my ($file_name,$dir,$suffix)=fileparse($file,qr/(\.fasta|\.fastq)\.*[^.]*/);
    
    if($no_R_match){ 
	if($file_name =~ /(.+)_([1|2])$/){
		$paired_files{$1.$suffix}[$2-1]=$file;
	} else {
		warn "Input file does not contain '_1'or '_2' in name: $file. Are you sure you should be using the --no_R_match option?";
	}

    } else {

    	if($file_name =~ /(.+)_R([1|2])[\.|_]*/){
		$paired_files{$1.$suffix}[$2-1]=$file;
        } else {
		warn "Input file does not contain '_R1_'and '_R2_' (or '_R1[.]' and '_R2[.]' in name: $file";
        }
   }
}

foreach my $name (keys %paired_files){

    # Unless full_name option set take name to be first field after delimiting by "_".
    my $out_name = undef;

    if(! $full_name){
  	my ($file_name,$dir,$suffix)=fileparse($name,qr/(\.fasta|\.fastq)\.*[^.]*/);
        $out_name = (split(/_/, $file_name))[0] . $suffix;
    } else {
	$out_name = $name;
    }

    my $pid = $pm->start and next; 

    my $cmd=join(' ','cat',@{$paired_files{$name}});
    
    my $out_file=$out_dir.'/'.$out_name;
    $cmd.= " > $out_file";
    print $cmd,"\n";
    unless($write_only){
	system($cmd);
    }
    $pm->finish;
}

__END__

=head1 Name

concat_paired_end.pl - Concatenate paired end files (fasta or fastq) (gzipped or not). By default the sample name is taken to be the first field after delmiting the file by "_".

=head1 USAGE

concat_paired_end.pl [-p [<# proc>] --write_only --no_R_match -h] -o <out_dir> <list of paired end files>

E.g.

#Note: Files must have "_R1_" and "_R2_" within the file name (or secondarily "_R1." and "_R2.")

concat_paired_end.pl sample1_R1_001.fastq sample1_R2_001.fastq sample2_R1_001.fastq sample2_R2_001.fastq

#Shorter way to do the same thing

concat_paired_end.pl *.fastq

#Specify alternate location for output files (instead of default current directory)

concat_paired_end.pl -o concatenated_reads *.fastq

#Run in parallel and use all CPUs

concat_paired_end.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

concat_paired_end.pl *.fastq -p 2

=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all  output files.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-w, --write_only>

Write the commands that will be run without actually running them.

=item B<-n, --no_R_match>

Identify forward and reverse pairs based on "_1" and "_2" directly before the suffix (e.g. _1.fastq) rather than "_R1_" and "_R2_".

=item B<-f, --full_name>

Sample names should be specified by full string upstream of forward/reverse identifier rather than just the first field when delimiting by "_".

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<concat_paired_end.pl> 

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

