#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;


my ($gzip,$write_only,$parallel,$help);
my $out_dir='./';
my $res = GetOptions("out_dir=s" => \$out_dir,
		     "parallel:i"=>\$parallel,
		     "gzip"=>\$gzip,
		     "write_only"=>\$write_only,
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
    if($file_name =~ /(.+)_R([1|2])[\.|_]*/){
	$paired_files{$1.$suffix}[$2-1]=$file;

### removed since fastqs often have "_1" and "_2" in ID 
###    #attempt different naming scheme
###    }elsif($file_name =~ /(.+)_([1|2])/){
####	$paired_files{$1.$suffix}[$2-1]=$file;

    }else{
###	warn "Input file does not contain '_R1_'and '_R2_' or '_1' and '_2' in name: $file";
	warn "Input file does not contain '_R1_'and '_R2_' (or '_R1[.]' and '_R2[.]' in name: $file";
    }
}

#my @out_files = map ($metaphlan_out_dir.$_, keys %paired_files);

foreach my $name (keys %paired_files){
    my $pid = $pm->start and next; 

    my $cmd=join(' ','cat',@{$paired_files{$name}});
    
    my $out_file=$out_dir.'/'.$name;
    $cmd.= " > $out_file";
    print $cmd,"\n";
    unless($write_only){
	system($cmd);
    }
    $pm->finish;
}

__END__

=head1 Name

concat_paired_end.pl - Concatenate paired end files (fasta or fastq) (gzipped or not)

=head1 USAGE

concat_paired_end.pl [-p [<# proc>] -h] -o <out_dir> <list of paired end files>

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

=item B<-h, --help>

Displays the entire help documentation.

=back

=head1 DESCRIPTION

B<concat_paired_end.pl> 

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

