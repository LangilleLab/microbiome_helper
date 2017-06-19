#!/usr/bin/env python

from __future__ import print_function
import argparse
import csv

__author__ = "Gavin Douglas"
__credits__ = ["Gavin Douglas"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Gavin Douglas"
__email__ = "gavin.douglas@dal.ca"
__status__ = "Development"

parser = argparse.ArgumentParser(description="Fix STAMP-formatted OTU table\
 so that all \"unidentified\" taxa are children of unique\
 parents. This is done by appending \"X\" to the name of the closest\
 higher-order taxonomic label that is defined. Additional \"X\" characters\
 will be added to distinguish undefined children at different levels.",

                                 epilog='''Usage example: fix_ITS2_spf.py -i\
otu_table.spf -o otu_table_fixed.spf''',

                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input", help="Input STAMP file", required=True)

parser.add_argument("-o", "--output", help="Output fixed STAMP file",
                    required=True)

parser.add_argument("-c", "--col_count", help="Number of taxonomic columns in\
                    file (default=7)", required=False, default=7, type=int)


def main():

    header_marker = 0

    args = parser.parse_args()

    outfile = open(args.output, "w")
    out_writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")

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

            # If no labels are unidentified then
            # print out line and move to next one.
            if "unidentified" not in line.lower():

                outfile.write(line + "\n")
                continue

            # Get list of all taxonomic levels.
            taxa = line_split[0:args.col_count]

            # Inititalize output list, string containing the most recently
            # defined parent, and a counter of the number of unidentified
            # levels since last defined parent.
            out_taxa = []
            defined_parent = "Unknown"
            num_unidentified = 0

            # Loop over each label (from higher to lower levels).
            # If label is unidentified then add Xs to end.
            # The number of Xs added is equal to the number of undefined
            # parents + 1.

            for label in taxa:
                if "unidentified" in label.lower():

                    num_unidentified += 1
                    out_taxa.append(defined_parent + "_" +
                                    num_unidentified * "X")
                else:
                    # If defined label then set as last defined parent and
                    # reset unidentified counter.
                    out_taxa.append(label)
                    defined_parent = label
                    num_unidentified = 0

            line_split[0:args.col_count] = out_taxa
            out_writer.writerow(line_split)

    outfile.close()


if __name__ == "__main__":
    main()
