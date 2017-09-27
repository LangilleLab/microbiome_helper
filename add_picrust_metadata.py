#!/usr/bin/python 

import sys
import argparse

parser = argparse.ArgumentParser(
                        description="Basic script to output function metadata "
                                    "in same order as predicted table. "
                                    "Metadata file should contain 1 function "
                                    "per line.",
                                    epilog="Usage example:\n"
                                    "add_picrust_metadata.py -i "
                                    "predicted.tab -m metadata.txt"
                                    "-o metadata_out.txt",
                        formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input", help="Predicted trait table.", 
                    required=True, type=str)
	
parser.add_argument("-m", "--meta", help="File containing metadata.", 
                    required=True, type=str)

parser.add_argument("-o", "--output", help="Output metadata line to concatenate.", 
                    required=True, type=str)

args = parser.parse_args()

count_file = args.input
metadata_file = args.meta
new_file = args.output

out_fh=open(new_file,'w')

def load_metadata(metadata_file_name):
    
    metadata={}
    metadata_fh=open(metadata_file_name)
    meta_name=metadata_fh.readline().strip().split('\t')[1]
    for line in metadata_fh:
        meta_fields=line.strip().split('\t')
        if len(meta_fields)> 1:
            metadata[meta_fields[0]]=meta_fields[1]
    return metadata,meta_name
        

count_fh=open(count_file)
header=count_fh.readline()
fields=header.strip().split('\t')
fields=fields[1:-1]

metadata,meta_name=load_metadata(metadata_file)

meta_line=[meta_name]
for field in fields:
    if field in metadata:
        if metadata[field] =='':
            meta_line.append('None')
        else:
            meta_line.append(metadata[field])
    
    else:
        print "Couldn't find metadata for "+field
        meta_line.append('None')

meta_line_str = "\t".join(meta_line)+"\n"
out_fh.write(meta_line_str)

