import argparse
import os
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(

            description="",

                epilog='''
                Examples for running:
                python filter_hits.py -hf hmmsearch_SCG_17_strains/SerS-Marinobacter_sp_v034.out -g SerS -sn Marinobacter_sp_v034 -o hmm_hits/SerS-Marinobacter_sp_v034.hits
                This would search the file hmmsearch_SCG_17_strains/SerS-Marinobacter_sp_v034.out for hits below the default e_value of e-25 and will write the number of hits to hmm_hits/SerS-Marinobacter_sp_v034.hits using the gene name SerS and sequence name Marinobacter_sp_v034

                python filter_hits.py -hf hmmsearch_SCG_17_strains/SerS-Marinobacter_sp_v034.out -g SerS -sn Marinobacter_sp_v034 -o hmm_hits/SerS-Marinobacter_sp_v034.hits -f --in_fasta genomes_17_strains/Marinobacter_sp_v034.ffn --out_fasta hmm_hits/SerS-Marinobacter_sp_v034.fasta
                This will do the same as above, but it will additionally find the sequences for the hits in the in_fasta file genomes_17_strains/Marinobacter_sp_v034.ffn and will write them to the out_fasta file hmm_hits/SerS-Marinobacter_sp_v034.fasta
                ''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-hf', '--hit_file', required=True,
                            type=str, help='output file from hmmsearch')
parser.add_argument('-g', '--gene_name', required=True,
                            type=str, help='Name of gene required for output')
parser.add_argument('-sn', '--sequence_name', required=True,
                            type=str, help='Name of sequence file required for output')
parser.add_argument('-o', '--out_file', required=True,
                            type=str, help='Name of output file, required for writing number of hits to')
parser.add_argument('-e', '--e_value', default='-25',
                            type=str, help='Threshold to include hits below this limit (default: %(default)s)')
parser.add_argument('-f', '--fasta', default=False, action='store_true',
                            help='Whether to write hits to a fasta file. If this is included, in_fasta and out_fasta must be set (default: %(default)s)')
parser.add_argument('--in_fasta', type=str, help='Input fasta file. Must be set if --fasta True')
parser.add_argument('--out_fasta', type=str, help='File name to write hmm hits to. Must be set if --fasta True')

args = parser.parse_args()
hit_file, e_value, fasta, in_fasta, out_fasta, gene_name, sequence_name, out_file = args.hit_file, args.e_value, args.fasta, args.in_fasta, args.out_fasta, args.gene_name, args.sequence_name, args.out_file

def check_file(fn):
    if not os.path.exists(fn):
        sys.exit("This file doesn't exist: "+fn)
    return

files_to_check = [hit_file]

# check whether in_fasta and out_fasta files have been given if fasta = True
if fasta:
    if in_fasta == None:
        sys.exit("You need to give a file for --in_fasta if fasta=True")
    if out_fasta == None:
        sys.exit("You need to give a name for --out_fasta if fasta=True")
    files_to_check.append(in_fasta)

if out_file == None:
    sys.exit("You need to give a name for --out_file")

# check files exist
for f in files_to_check:
    check_file(f)

# get hits
e_value = 10**(int(e_value))
hits = []
for row in open(hit_file, 'r'):
    if row[0] == '#': continue
    row = list(filter(None, row.replace('\n', '').split(' ')))
    if float(row[4]) <= e_value:
        hits.append(row[0])

# write number of hits to out_file, using sequence_name and gene_name as two other columns in this output
with open(out_file, 'w') as f:
    w = f.write(gene_name+'\t'+sequence_name+'\t'+str(len(hits))+'\n')

# get the hits from the input fasta and write these to the output fasta
if fasta:
    hit_seqs = []
    for record in SeqIO.parse(in_fasta, "fasta"):
        if record.id in hits:
            hit_seqs.append(record)
    
    w = SeqIO.write(hit_seqs, out_fasta, "fasta")
