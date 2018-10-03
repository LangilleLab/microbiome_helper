#!/usr/bin/env python

from __future__ import print_function, division
from tempfile import TemporaryDirectory
from multiprocessing import Pool
import argparse
import pandas as pd
import subprocess
import os

__author__ = "Gavin Douglas"

parser = argparse.ArgumentParser(
                        description="Creates table of number of reads per " + \
                                    "sample from QIIME2 qza files. Column " + \
                                    "names will be QZA filenames with " + \
                                    "extension removed.",
                        formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('QZA_files', metavar='QZA', type=str, nargs='+',
                    help='QIIME2 artifacts containing gzipped fastq files.')


parser.add_argument("-p", "--proc", required=False, type=int, default=1,
                    help="Number of processes to run in parallel") 

parser.add_argument("-o", "--output", help="Name of output read count table",
                    required=True, type=str)

def fastq_read_num(fastq):
    '''Runs system calls to determine number of reads in gzipped fastq''' 

    zcat_process = subprocess.Popen(["zcat", fastq], stdout=subprocess.PIPE)

    wc_process = subprocess.Popen(["wc", "-l"], stdin=zcat_process.stdout,
                                  stdout=subprocess.PIPE)

    zcat_process.stdout.close()

    wc_output, wc_err = wc_process.communicate()

    return(int(int(wc_output.decode('ascii'))/4))


def process_fastq(fastq):
    '''For an individual FASTQ file determine the sample name, whether it's R1
    or R2, and the number of reads in the file. Returns list of the sample and
    orientation joined by "_" as first element and the number of reads as
    second element.'''

    fastq_basename = os.path.basename(fastq)

    sample = fastq_basename.split("_")[0]

    if "_R1_" in fastq_basename:
        read_id = "R1"
    elif "_R2_" in fastq_basename:
        read_id = "R2"
    else:
        raise ValueError("Neither _R1_ or _R2_ in filename " + fastq_basename)

    return(["_".join([sample,read_id]),
            fastq_read_num(fastq)])


def reads_per_qza_gzipped_fastq(qza, num_proc):
    '''Decompress QZA of FASTQs and return dictionary of read counts per
    sample for both R1 and R2 fastqs (R2 only if applicable).'''

    sample_reads = {}

    # In temporary directory decompress the qza and output all FASTQs.
    with TemporaryDirectory() as tmp:

        qza_name = os.path.basename(qza)
        sample_reads[qza_name] = {}

        export_cmd = ["qiime", "tools", "export", "--input-path", qza, "--output-path", tmp]
        subprocess.call(export_cmd)

        # Get list of all FASTQs in tmp folder.
        tmp_fastqs = [os.path.join(tmp, f) for f in os.listdir(tmp) if os.path.isfile(os.path.join(tmp, f)) and ".fastq.gz" in f]

        pool_set = Pool(num_proc)

        sample_reads_raw = pool_set.map(process_fastq, tmp_fastqs)

        # Add these output values to a dictionary to be returned.
        for raw_out in sample_reads_raw:
            sample_reads[qza_name][raw_out[0]] = raw_out[1]

    return(sample_reads)


def main():

    args = parser.parse_args()

    # Initialize list containing dictionaries for each QZA.
    sample_reads_dicts = []

    # Loop over all input QZA files.
    for qza in args.QZA_files:
        sample_reads_dicts += [reads_per_qza_gzipped_fastq(qza, args.proc)]

    # Combine all dictionaries into one.
    sample_reads_combined_dict = {}

    for tab_dict in sample_reads_dicts:
        sample_reads_combined_dict.update(tab_dict)

    # Make pandas dataframe from thsi combined dictionary.
    sample_reads = pd.DataFrame.from_dict(sample_reads_combined_dict, orient='columns')

    # Re-order indices and columns alphabetically.
    sample_reads.sort_index(axis=0, inplace=True)
    sample_reads.sort_index(axis=1, inplace=True)

    sample_reads.to_csv(path_or_buf=args.output,  sep="\t",
                        index_label="dataset", na_rep="NaN")

if __name__ == "__main__":
    main()
