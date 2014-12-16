#!/usr/bin/env python
from __future__ import division

__author__ = "Morgan Langille"
__credits__ = ["Morgan Langille"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Morgan Langille"
__email__ = "morgan.g.i.langille@gmail.com"
__status__ = "Development"


from cogent.util.option_parsing import parse_command_line_parameters, make_option
from biom.parse import parse_biom_table
from os.path import join,splitext
import gzip
import sys
script_info = {}
script_info['brief_description'] = "Convert a PICRUST KO BIOM file to a compatible STAMP profile table FOR A GIVEN KEGG PATHWAY"
script_info['script_description'] = "This handles the one-to-many mapping of KOs to KEGG Pathways. Since STAMP does not handle one-to-many mappings, the -p option should be used to limit to only a single KEGG Pathway. Also, this pulls the KEGG Ortholog descriptions for each KO id."

script_info['script_usage'] = [\
("Minimum Requirements (not compatible with STAMP)","","%prog -i ko.biom  > ko.spf"),
("Using PICRUSt KO table","","%prog -i ko.biom -p 'Benzoate_degradation' > ko_benzoate.spf"),
]

script_info['output_description']= "A STAMP profile containing KOs for a given KEGG Pathway"

script_info['optional_options']=[\
make_option('-p','--pathway',default=None,type="string",help='Name of KEGG Pathway.[default: %default')]

script_info['required_options'] = [\
make_option('-i','--input_fp',type="existing_filepath",help='the input PICRUSt KO filepath in biom format (.biom)')]

script_info['version'] = __version__
       

def main():
    option_parser, opts, args =\
                   parse_command_line_parameters(**script_info)

    file_name=opts.input_fp
    metadata_name='KEGG_Pathways'
    id_description='KEGG_Description'

    #allow file to be optionally gzipped (must use extension '.gz')
    ext=splitext(file_name)[1]
    if (ext == '.gz'):
        table = parse_biom_table(gzip.open(file_name,'rb'))
    else:
        table = parse_biom_table(open(file_name,'U'))

    #make the header line
    header=['Level_1','Level_2','KEGG_Pathway','KEGG_Ortholog']
    
    #add the sample ids to the header line
    header.extend(table.SampleIds)

    print "\t".join(header)
    
    max_len_metadata=3
    #now process each observation (row in the table)
    for obs_vals,obs_id,obs_metadata in table.iterObservations():
        #import pdb;pdb.set_trace()
        ko_description=obs_metadata[id_description]
        if type(ko_description) is list:
            ko_description="; ".join(ko_description)

        for pathway in obs_metadata[metadata_name]:
           
            #Check to make sure the metadata has the correct number of levels
            
            if len(pathway) < max_len_metadata:
                for i in range(max_len_metadata - len(pathway)):
                    pathway.append('')

            #skip pathways that don't match what the user wants to keep
            if opts.pathway:
                if not opts.pathway == pathway[2]:
                    continue

            ko_label=obs_id+': '+ko_description
            
            #Add the ko_label as the last "Level"
            pathway.append(ko_label)

            #Add count data to the row
            pathway.extend(map(str,obs_vals))
            print "\t".join(pathway)
        
    

if __name__ == "__main__":
    main()

