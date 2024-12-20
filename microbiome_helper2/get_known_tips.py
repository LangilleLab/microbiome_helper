from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(

            description="",

                epilog='''
                    Add examples here
                    ''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--in_fasta', required=True,
                            type=str, help='Input fasta to take tip names from')
parser.add_argument('-o', '--out_file', required=True,
                            type=str, help='Output file name to save csv of sequence names to')

args = parser.parse_args()

in_fasta, out_file = args.in_fasta, args.out_file

names = []
for record in SeqIO.parse(in_fasta, "fasta"):
    names.append(record.id)

with open(out_file, "w") as f:
    for name in names:
        w = f.write(name+'\n')
