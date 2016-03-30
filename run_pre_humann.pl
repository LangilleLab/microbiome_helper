#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;

my ($parallel,$help);
my $out_dir='./';
my $db;
my $cluster;
my $tmp_dir;
my $write_only;
my %search_types=('diamond'=>1,'blast'=>1);
my $search_type='diamond';
my $diamond_location='diamond';
my $force;
my $res = GetOptions("out_dir=s" => \$out_dir,
		     "parallel:i"=>\$parallel,
		     "db=s" =>\$db,
		     "search_type=s"=>\$search_type,
		     "cluster_submission"=>\$cluster,
		     "tmp_dir=s"=>\$tmp_dir,
		     "location=s"=>\$diamond_location,
		     "write_only"=>\$write_only,
		     "force"=>\$force,
		     "help"=>\$help,
    )or pod2usage(2);

pod2usage(-verbose=>2) if $help;

my @files=@ARGV;

pod2usage($0.': You must provide a list of fasta/fastq files to be searched against KEGG.') unless @files;

pod2usage("Your --search_type $search_type is not one ",join(keys %search_types)) unless $search_types{$search_type}; 

#set default for tmp dir (unless already set by the user)
unless($tmp_dir){
    #if submitting jobs to computer cluster
    if($cluster){
	#use env variable $TMPDIR
	$tmp_dir='\$TMPDIR';
    }else{
	$tmp_dir='/tmp';
    }
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
foreach my $file (@files){
    my ($file_name,$dir)=fileparse($file, qr/\.[^.]*/);

    my $out_file=$out_dir.'/'.$file_name;
    my $out_file_daa=$out_file.'.daa';
    my $out_file_txt=$out_file.'.txt';

    my $cmd;
    if($search_type eq 'diamond'){
	$db="/home/shared/kegg/diamond_db/kegg.reduced" unless $db;
	$cmd="$diamond_location blastx -p $cpu_count -d $db -q $file -a $out_file -t $tmp_dir -c 1;$diamond_location view -a $out_file_daa -o $out_file_txt;rm $out_file_daa";
    }elsif($search_type eq 'blast'){
	$db="/home/shared/kegg/blast_db/kegg.reduced" unless $db;
	$cmd="blastx -num_threads $cpu_count -outfmt 6 -db $db -query $file -out $out_file_txt"
    }

    if($cluster){
	my $log_dir='./log/';
	`mkdir -p $log_dir`;
	my $stderr=$log_dir.$file_name.'.stderr';
	my $stdout=$log_dir.$file_name.'.stdout';
	$cmd = "echo \"echo \\\$TMPDIR; hostname; time $cmd\" | qsub -l h_rt=46:00:00 -pe openmp $cpu_count -l h_vmem=4G -V -cwd -e \"$stderr\" -o \"$stdout\" -N \"$file_name\" -S /bin/bash";
    }

    #only do run if output doesn't exist (or force option is used)
    if(!(-e $out_file_txt) || $force){
	print $cmd,"\n";
	unless($write_only){
	    system($cmd);
	}
    }
}

__END__

=head1 Name

run_pre_humann.pl - A simple wrapper for pre_humann to screen for human sequences in metagenomic data

=head1 USAGE

run_pre_humann.pl [-d <db> -p [<# proc>] -o <out_dir> -h] <list of fastq or fasta files>

E.g.

run_pre_humann.pl sample1_assembled.fastq sample2_assembled.fastq

#Shorter way to do the same thing

run_pre_humann.pl *.fastq

#Specify alternate location for output files (instead of default current directory)

run_pre_humann.pl -o screened_reads *.fastq

#Run in parallel and use all CPUs

run_pre_humann.pl *.fastq -p

#Run in parallel limit to only 2 CPUs

run_pre_humann.pl *.fastq -p 2


=head1 OPTIONS

=over 4

=item B<-o, --out_dir <file>>

The name of the output directory to place all output files.

=item B<-d, --db <file>>

Path to KEGG database. 
Default if running DIAMOND: "/home/shared/kegg/diamond_db/kegg.reduced".
Default if running BLAST: "/home/shared/kegg/blast_db/kegg.reduced".

Note that these databases may need to be generated with same version of DIAMOND/BLAST that you are using. 

=item B<-s, --search_type <blast|diamond>>

Specifies either "blast" or "diamond" (default: "diamond").

=item B<-c, --cluster_submission>

Set flag if running on cluster (when set == TRUE). 

=item B<-t, --tmp_dir <path>>

Path to tmp directory to use. By default is /tmp except on clusters where the TMP environment variable is used by default.

=item B<-l, --location <path>>

Path to DIAMOND binary. By default is "diamond", so it assumed to be in your $PATH. 

=item B<-w, --write_only>

Flag that will result in this script printing the commands to be run, but wont actually run them.

=item B<-f, --force>

Flag to indicate that old output files should be overwritten.

=item B<-p, --parallel [<# of proc>]>

Using this option without a value will use all CPUs on machine, while giving it a value will limit to that many CPUs. Without option only one CPU is used. 

=item B<-h, --help>

Displays the entire help documentation.

=back

		     
=head1 DESCRIPTION

B<run_pre_humann.pl> This script allows for more automated running of the pre_humann program on multiple fastq files. 
 
Before use make sure you have installed the "pre_humann.pl" (along with corresponding databases) program so it is accesible from your PATH

=head1 AUTHOR

Morgan Langille, E<lt>morgan.g.i.langille@gmail.comE<gt>

=cut

