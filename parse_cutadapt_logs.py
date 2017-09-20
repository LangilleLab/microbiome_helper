#!/usr/bin/env python

from __future__ import print_function

__author__ = "Gavin Douglas"
__credits__ = ["Gavin Douglas"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Gavin Douglas"
__email__ = "gavin.douglas@dal.ca"
__status__ = "Development"

import argparse
import re
import os

parser = argparse.ArgumentParser(description="Parse cutadapt logfiles to \
produce a table of read counts per sample.", epilog='''Usage example:

parse_cutadapt_logs.py -i primer_trimmed_fastqs -m "_log.txt" \
-o cutadapt_log.txt

''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input", help="Folder containing cutadapt \
logfiles", required=True, type=str)

parser.add_argument("-m", "--match", help="String to match cutadapt log \
filenames", required=False, default="_log.txt", type=str)

parser.add_argument("-s", "--sample_delim", help="String to split filename on \
to determine sample name (assumed to be 1st field after splitting",
                    required=False, default="_", type=str)

parser.add_argument("-o", "--output", help="Name of output read count table",
                    required=False, default="cutadapt_log.txt", type=str)


def main():

    args = parser.parse_args()

    # Check whether input folder exists.
    if not os.path.isdir(args.input):
        raise Exception("Input path " + args.input + " is not a directory")

    # Read in cutadapt log filenames.
    logfiles = [f for f in os.listdir(args.input) if args.match in f]

    # Throw exception if no logfiles matched string.
    if len(logfiles) == 0:
        raise Exception("No logfiles found in folder " +
                        args.input + " matching string " + args.match)

    # Open output file for writing.
    outfile = open(args.output, "w")

    # Put samples in set to double-check they are all unique.
    past_samples = set()

    # Set flag for whether header printed or not.
    header_printed = False

    # Loop over all logfiles.
    for logfile in logfiles:

        sample = logfile.split(args.sample_delim)[0]

        # Check whether sample already seen and if not add to set.
        if sample in past_samples:
            raise Exception("Parsed sample name " + sample +
                            " is not unique!")

        past_samples.add(sample)

        log = args.input + "/" + logfile

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
                for i in range(len(strings2match)):

                    if strings2match[i] in line:

                        count_regex = r".*" + re.escape(strings2match[i]) + \
                                       "\s+" + "([0-9]+,*[0-9]*).*"

                        count_match = re.match(count_regex, line).group(1)

                        if(count_match):
                            category_match = count_categories[i]

                            # Check whether this string was matched already.
                            if category_match in read_counts:
                                raise Exception("Multiple matches to string " +
                                                strings2match[i] + " in " +
                                                log)

                            # Remove comma.
                            count_match = count_match.replace(",", "")

                            read_counts[category_match] = count_match

        # If library was not defined then throw error.
        if not library:
            raise Exception("Sequencing library type could not be " +
                            "determined for file " + log)

        # Initialize list to write out per sample.
        out_counts = [sample]

        # Loop over strings/categories indices again and make sure they are
        # all present.
        for i in range(len(strings2match)):

            if count_categories[i] not in read_counts:
                raise Exception("Could not find line matching " +
                                strings2match[i] + " for count of " +
                                count_categories[i] + " in " + log)

            out_counts = out_counts + [read_counts[count_categories[i]]]

        print("\t".join(out_counts), file=outfile)

    outfile.close()


if __name__ == "__main__":
    main()
