#!/usr/bin/env python

from __future__ import print_function

__author__ = "Gavin Douglas"
__credits__ = ["Gavin Douglas"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Gavin Douglas"
__email__ = "gavin.douglas@dal.ca"
__status__ = "Development"

import argparse
import re
import os
import sys

parser = argparse.ArgumentParser(
                        description="Parse cutadapt logfiles to: produce a "
                                    "table of read counts per sample.",
                                    epilog="Usage example:\n"
                                    "parse_cutadapt_logs.py -i "
                                    "primer_trimmed_fastqs/*cutadapt_log.txt "
                                    "-o cutadapt_log.txt",
                        formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input", help="Cutadapt logfiles", required=True,
                    type=str, nargs='+')

parser.add_argument("-s", "--sample_delim",
                    help="String to split filename on to determine sample "
                         "name (assumed to be 1st field after splitting)",
                    required=False, default="_", type=str)

parser.add_argument("-o", "--output", help="Name of output read count table",
                    required=False, default="cutadapt_log.txt", type=str)


def main():

    args = parser.parse_args()

    # Open output file for writing.
    outfile = open(args.output, "w")

    # Put samples in set to double-check they are all unique.
    past_samples = set()

    # Set flag for whether header printed or not.
    header_printed = False

    # Loop over all logfiles.
    for log in args.input:

        logfile = os.path.basename(log)

        sample = logfile.split(args.sample_delim)[0]

        # Check whether sample already seen and if not add to set.
        if sample in past_samples:
            sys.exit("Parsed sample name {sample} is not unique!".format(
                                                                sample=sample))

        past_samples.add(sample)

        # Initialize library object which will be either "single" or "paired"
        library = None

        # Initialize empty dict which will be actual counts parsed.
        read_counts = {}

        with open(log, "r") as infile:

            for line in infile:

                # Check if library defined. If not then check if library info
                # is on line. If not then skip line. Set strings that should
                # be matched as well depending on library type and what
                # categories these strings correspond to.
                if not library:
                    if "in paired-end mode ..." in line:

                        library = "paired"

                        strings2match = ["Total read pairs processed:",
                                         "Read 1 with adapter:",
                                         "Read 2 with adapter:",
                                         "Pairs written (passing filters):"]

                        count_categories = ["pairs_in", "read1_match",
                                            "read2_match", "pairs_passed"]

                    elif "in single-end mode ..." in line:

                        library = "single"

                        strings2match = ["Total reads processed:",
                                         "Reads with adapters:",
                                         "Reads written (passing filters):"]

                        count_categories = ["reads_in", "reads_match",
                                            "reads_passed"]

                    # Check if header printed yet and if not then do so.
                    if library and not header_printed:
                        header_out = ["sample"] + count_categories
                        print("\t".join(header_out), file=outfile)
                        header_printed = True

                    # Skip to next line.
                    continue

                # If library defined then check for matching strings.
                for string2match, count_category in zip(strings2match,
                                                        count_categories):

                    if string2match in line:

                        count_regex = r".*" + re.escape(string2match) + \
                                       "\s+" + "([0-9,]+).*"

                        count_full_match = re.match(count_regex, line)

                        if count_full_match:

                            count_match = count_full_match.group(1)

                            # Check whether this string was matched already.
                            if count_category in read_counts:
                                sys.exit("Multiple matches to string "
                                         "{string2match} in {log}".
                                         format(string2match=string2match,
                                                log=log))

                            # Remove comma if present.
                            count_match = count_match.replace(",", "")

                            read_counts[count_category] = count_match

        # If library was not defined then throw error.
        if not library:
            sys.exit("Sequencing library type could not be "
                     "determined for file {log}".format(log=log))

        # Initialize list to write out per sample.
        out_counts = [sample]

        # Loop over strings/categories indices again and make sure they are
        # all present.
        for string2match, count_category in zip(strings2match,
                                                count_categories):

            if count_category not in read_counts:
                sys.exit("Could not find line matching "
                         "{string2match} for count of "
                         "{count_category} in {log}".format(
                                                string2match=string2match,
                                                count_category=count_category,
                                                log=log))

            out_counts = out_counts + [read_counts[count_category]]

        print("\t".join(out_counts), file=outfile)

    outfile.close()


if __name__ == "__main__":
    main()
