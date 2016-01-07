#!/usr/bin/env python
from __future__ import division

__author__ = "Morgan Langille"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"

import argparse
from os.path import join
import sys
import re
#Requires BIOM v2.1
from biom import load_table
from biom.util import biom_open,HAVE_H5PY
import math

parser = argparse.ArgumentParser(description="Remove OTUs with low confidence, based on known 0.1% bleed through between runs on the Illumina MiSeq. Cufoff=(# of seqs)/(# of samples)*0.001.", 
                                 epilog='''Examples of Usage:
#OTU table from QIIME:
remove_low_confidence_otus.py -i otu_table.biom -o otu_table_clean.biom

'''
                                 ,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input",help="Input BIOM file", required=True)

parser.add_argument("-o","--output",help="Output BIOM file", required=True)


def write_biom_table(biom_table, biom_table_fp, compress=True,
                     write_hdf5=HAVE_H5PY, format_fs=None):
    """Writes a BIOM table to the specified filepath

    Parameters
    ----------
    biom_table : biom.Table
        The table object to write out
    biom_table_fp : str
        The path to the output file
    compress : bool, optional
        Defaults to ``True``. If True, built-in compression on the output HDF5
        file will be enabled. This option is only relevant if ``write_hdf5`` is
        ``True``.
    write_hdf5 : bool, optional
        Defaults to ``True`` if H5PY is installed and to ``False`` if H5PY is
        not installed. If ``True`` the output biom table will be written as an
        HDF5 binary file, otherwise it will be a JSON string.
    format_fs : dict, optional
        Formatting functions to be passed to `Table.to_hdf5`

    Notes
    -----
    This code was adapted from QIIME 1.9
    """
    generated_by = "Microbiome Helper"

    if write_hdf5:
        with biom_open(biom_table_fp, 'w') as biom_file:
            biom_table.to_hdf5(biom_file, generated_by, compress,
                               format_fs=format_fs)
    else:
        with open(biom_table_fp, 'w') as biom_file:
            biom_table.to_json(generated_by, biom_file)


def main():
    args = parser.parse_args()

    input_filename = args.input
    table = load_table(input_filename)

    output_filename=args.output

    num_samples=table.length(axis='sample')
    num_reads=table.sum(axis='whole')
    error_rate=0.001

    cut_off=math.ceil(num_reads/num_samples*error_rate)

    def filter_fn(val, id_, md):
        return (val.sum() > cut_off)

    table.filter(filter_fn, axis='observation', inplace=True)

    write_biom_table(table, output_filename)

    print "OTUs with more than %d sequences were kept." % (cut_off)
    

if __name__ == "__main__":
    main()

