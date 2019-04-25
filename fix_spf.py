#!/usr/bin/env python

from __future__ import print_function
import argparse
import csv
import re

__author__ = "Gavin Douglas"
__license__ = "GPL"
__version__ = "0.1"

parser = argparse.ArgumentParser(

description="Fix STAMP-formatted OTU table outputted by QIIME2 with SILVA "
            "taxonomy. Will convert all taxa containing \"uncultured\", "
            "\"Ambiguous_taxa\", \"metagenome\", or \"unidentified\", with "
            "\"Unclassified\". Will then replace all taxa containing "
            "\"unknown\" with the preceeding taxonomic label followed by "
            "\"X\".",

epilog='''Usage example: fix_spf.py -i otu_table.spf -o otu_table_fixed.spf''',

formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input", help="Input STAMP file", required=True)

parser.add_argument("-o", "--output", help="Output fixed STAMP file",
                    required=True)

parser.add_argument("-c", "--col_count", required=False, default=7, type=int,
                    help="Number of taxonomic columns in file (default=7)", )


def main():

    header_marker = 0

    args = parser.parse_args()

    outfile = open(args.output, "w")
    out_writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")

    str2replace = ['uncultured', 'ambiguous_taxa', 'metagenome',
                   'unidentified']

    with open(args.input, "r") as infile:

        for line in infile:

            # Remove newline character and split by tab.
            line = line.rstrip("\r\n")
            line_split = line.split("\t")

            # Print out header and go to next line.
            if header_marker == 0:
                header_marker += 1
                outfile.write(line + "\n")
                continue

            # If no labels are in the set to replace then print out line and
            # move to next one.
            str_matches = 0

            for s in str2replace + ['unknown']:
                if s in line.lower():
                    str_matches += 1

            if str_matches == 0:
                outfile.write(line + "\n")
                continue

            # Get list of all taxonomic levels.
            taxa = line_split[0:args.col_count]

            # Loop through taxa and replace any ids in set of strings to
            # replace with "Unclassified". For any taxa containing "unknown",
            # replace these ids with the preceeding taxonomic level, but with
            # the correct DX level and followed by "X".
            out_taxa = []

            # Loop over each label (from higher to lower levels).
            for label_i, label in enumerate(taxa):
                str_match = False

                for s in str2replace:
                    if s in label.lower():
                        str_match = True

                if str_match:
                    out_taxa.append('Unclassified')
                
                elif 'unknown' in label.lower():
                    pre_label_i = label_i - 1
                    if pre_label_i < 0 or taxa[pre_label_i] == 'Unclassified':
                        out_taxa.append('Unclassified')
                    else:
                        pre_label_search = re.search("D_\d+__(.*)", taxa[pre_label_i])
                        pre_label_taxon = pre_label_search.group(1)
                        label_level = re.match("D_\d+__", label).group(0)
                        out_taxa.append(label_level + pre_label_taxon + "X")

                else:
                    out_taxa.append(label)

            line_split[0:args.col_count] = out_taxa
            out_writer.writerow(line_split)

    outfile.close()


if __name__ == "__main__":
    main()
