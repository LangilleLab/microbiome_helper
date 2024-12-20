from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser(

            description="",

                epilog='''
                    Add examples here
                        check_alignment.py
                        ''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--ref_msa', required=True,
                            type=str, help='Reference multiple sequence alignment')
parser.add_argument('--in_msa', required=True,
                            type=str, help='Input multiple sequence alignment - it is expected that this will contain query sequences as well as reference sequences')
parser.add_argument('--in_fasta', required=True,
                            type=str, help='Unaligned fasta file containing only the query sequences')
parser.add_argument('--out_ref_msa', required=True,
                            type=str, help='Output reference multiple sequence alignment')
parser.add_argument('--out_msa', required=True,
                            type=str, help='Output multiple sequence alignment containing filtered query sequences')

args = parser.parse_args()

ref_msa, in_msa, in_fasta, out_ref_msa, out_msa = args.ref_msa, args.in_msa, args.in_fasta, args.out_ref_msa, args.out_msa

def read_fasta(fn):
    records = {}
    for record in SeqIO.parse(fn, "fasta"):
        records[record.id] = str(record.seq)
    return(records)

def check_alignments(raw_seqs, aligned_seqs, min_align):
    poorly_aligned = []
    passing = {}
    min_study_seq_length = None
    max_study_seq_length = None
    raw_seqnames = set(raw_seqs.keys())
    for seq_id in raw_seqnames:
        orig_seq = raw_seqs[seq_id]
        aligned_seq = aligned_seqs[seq_id]
        orig_seq = orig_seq.replace("-", "")
        orig_seq = orig_seq.replace(".", "")
        aligned_seq = aligned_seq.replace("-", "")
        aligned_seq = aligned_seq.replace(".", "")
        orig_seq_length = len(orig_seq)
        if len(aligned_seq) < orig_seq_length * min_align:
            poorly_aligned.append(seq_id)
        else:
            passing[seq_id] = aligned_seqs[seq_id]
        if not min_study_seq_length or orig_seq_length < min_study_seq_length:
            min_study_seq_length = orig_seq_length
        if not max_study_seq_length or orig_seq_length > max_study_seq_length: 
            max_study_seq_length = orig_seq_length
    if len(poorly_aligned) == len(raw_seqnames):
        print("All sequences aligned poorly to reference sequences")
    elif len(poorly_aligned) > 0:
        print(str(len(poorly_aligned))+" of "+str(len(raw_seqnames))+" aligned poorly to reference sequences")
    if min_study_seq_length == max_study_seq_length:
        print("All raw input sequences were the same length")
    else:
        print("Raw input sequences ranged in length from "+str(min_study_seq_length)+" to "+str(max_study_seq_length))
    return(passing)

ref_seqnames = set(list(read_fasta(ref_msa).keys()))
study_seqs = read_fasta(in_fasta)
align_out = read_fasta(in_msa)

ref_align_subset = {seq: align_out[seq] for seq in ref_seqnames}
study_align_subset = check_alignments(raw_seqs=study_seqs, aligned_seqs=align_out, min_align=0.8)

ref_align = []
for seq in ref_align_subset:
    ref_align.append(SeqRecord(Seq(ref_align_subset[seq]), id=seq, description=''))

study_align = []
for seq in study_align_subset:
    study_align.append(SeqRecord(Seq(study_align_subset[seq]), id=seq, description=''))

SeqIO.write(ref_align, out_ref_msa, "fasta")
SeqIO.write(study_align, out_msa, "fasta")
