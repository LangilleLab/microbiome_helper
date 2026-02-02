#!/usr/bin/env python

import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(

    description="Add function descriptions to HUMAnN2-formatted function "
                "abundance table, in either stratified or unstratified format. "
                "Note that this script assumes that function descriptions are "
                "not already present. Input mapfile needs to be a tab-delimited "
                "file without a header and with two columns. The first column "
                "should contain function ids and the second column should "
                "contain descriptions. One pair of ids and descriptions per "
                "line.",
    epilog='''

Usage:
descrip_humann2.py -i table.tsv -m /path/to/pathway_descrip.tsv -o table_w_descrip.tsv --add new

''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', metavar='PATH', required=True, type=str,
                    help='Input function abundance table.')

parser.add_argument('-o', '--output', metavar='PATH', type=str, required=True,
                    help='Output function abundance table with added '
                         'description column. If the extension \".gz\" '
                         'is added the table will automatically be gzipped.')

parser.add_argument('-m', '--map', metavar='PATH', type=str,
                    help='An input map table linking function ids to '
                         'descriptions for each function. ')

parser.add_argument('-a', '--add', metavar='new|existing', type=str,
                    default='new', choices=['new', 'existing'],
                    help='Parameter to specify that description should be '
                         'added to existing function column or whether a new '
                         'column for the description (and the taxonomy info '
                         'if present) should be created. One of either '
                         '\"existing\" or \"new\" should be input. A new '
                         'column will be created by default.')

parser.add_argument('--missing_descrip', metavar='STRING', type=str,
                    default="NA",
                    help='String to fill in for cases where no description is '
                         'found for a function id.')

parser.add_argument('--missing_taxon', metavar='STRING', type=str,
                    default="NA",
                    help='String to fill in for stratified tables when no '
                         'taxon string is present.')

def main():

    args = parser.parse_args()
    function_tab = pd.read_csv(args.input, sep="\t", low_memory=False,
                               dtype={'function': str})

    map_tab = pd.read_csv(args.map, sep="\t", index_col=0, header=None,
                          names=["function", "description"],
                          low_memory=False, dtype=object)

    # Assume first column contains function ids (and potentially taxon links 
    # too)
    first_col = function_tab.iloc[: , 0]

    # Check to see if file is stratified or not.
    strat_check = first_col.str.contains(pat="\|").any()

    # If stratified table then split function and taxonomy.
    if strat_check:
        function_info = pd.DataFrame(columns=['function', 'description', 'taxon'])

        function_info[['function', 'taxon']] = function_tab['function'].str.split('|', expand=True)

        function_info.loc[function_info['taxon'].isnull(), 'taxon'] = args.missing_taxon

    else:
        function_info = pd.DataFrame(columns=['function', 'description'])

        function_info["function"] = first_col

    # Reindex mapfile to match order of functions in function table (and for
    # ids to be duplicated if the table is stratified.
    map_tab = map_tab.reindex(function_info.function,
                              fill_value=args.missing_descrip)

    function_info["description"] = list(map_tab["description"])

    if args.add == "existing":
        if strat_check:
            function_tab.iloc[: , 0] = function_info['function'] + ':' + function_info['description'] + '|' + function_info['taxon']
        else:
            function_tab.iloc[: , 0] = function_info['function'] + ':' + function_info['description']
    elif args.add == "new":
        function_tab.drop(function_tab.columns[0], axis=1, inplace=True)

        function_tab = pd.concat([function_info, function_tab], sort=False, axis=1)

    function_tab.to_csv(path_or_buf=args.output, sep="\t",
                        index=False, compression="infer")

if __name__ == "__main__":
    main()
