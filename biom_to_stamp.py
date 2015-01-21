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
import re

script_info = {}
script_info['brief_description'] = "Convert a BIOM table to a compatible STAMP profile table."
script_info['script_description'] = "Metadata will be parsed and used as hiearachal data for STAMP."

script_info['script_usage'] = [\
("Minimum Requirments","","%prog table1.biom > table1.spf"),
("OTU table from QIIME","","%prog -m taxonomy otu_table.biom > otu_table.spf"),
("KO file from PICRUSt","","%prog -m KEGG_Description ko.biom > ko.spf"),
("KEGG Pathways table from PICRUSt","","%prog -m KEGG_Pathways ko_L3.biom > ko_L3.spf"),
("Function table from MG-RAST","","%prog -m ontology table1.biom > table1.spf")
]

script_info['output_description']= "Output is written to STDOUT"

script_info['optional_options'] = [\
    make_option('-m','--metadata',default=None,type="string",help='Name of metadata. [default: %default]')]


script_info['disallow_positional_arguments'] = False

script_info['version'] = __version__
       
def process_metadata(metadata,metadata_name,obs_id):
    if metadata_name =='taxonomy':
        fixed_metadata=[]
        for idx,val in enumerate(metadata):
            if(re.match(r'[a-z]__$',val)):
                fixed_metadata.append("Unclassified")
            else:
                fixed_metadata.append(val)
        return fixed_metadata

    elif metadata_name == 'KEGG_Pathways':
        if metadata[0]=='Unclassified':
            #Remove "Unclassified" from the first of the levels
            del metadata[0]
            metadata.append(metadata[-1]+'_Unclassified')
        return metadata
    elif metadata_name == 'KEGG_Description':
        single_metadata= ' or '.join(metadata)
        single_metadata=obs_id+': '+single_metadata
        return [single_metadata]
    else:
        return metadata

    

    

def main():
    option_parser, opts, args =\
                   parse_command_line_parameters(**script_info)

    min_args = 1
    if len(args) < min_args:
       option_parser.error('A BIOM file must be provided.')

    file_name = args[0]

    #allow file to be optionally gzipped (must use extension '.gz')
    ext=splitext(file_name)[1]
    if (ext == '.gz'):
        table = parse_biom_table(gzip.open(file_name,'rb'))
    else:
        table = parse_biom_table(open(file_name,'U'))

    metadata_name=opts.metadata

    if metadata_name is None:
        max_len_metadata=0
    elif metadata_name == 'KEGG_Description':
        max_len_metadata=1
    elif table.ObservationMetadata and metadata_name in table.ObservationMetadata[0]:
       
        max_len_metadata = max(len(p[metadata_name]) for p in table.ObservationMetadata)
    else:
        raise ValueError("'"+metadata_name+"' was not found in the BIOM table. Please try changing --metadata to a valid metadata field.")

    include_obs_id=True
    if metadata_name in ["KEGG_Pathways","KEGG_Description",'taxonomy']:
        include_obs_id=False

    #make the header line
    header=[]

    #make simple labels for each level in the metadata (e.g. 'Level_1', 'Level_2', etc.) "+1" for the observation id as well.
    extra_levels=0
    if include_obs_id:
        extra_levels=1
        
    for i in range(max_len_metadata+extra_levels):
        header.append('Level_'+ str(i+1))
    
    #add the sample ids to the header line
    header.extend(table.SampleIds)
    
    print "\t".join(header)

    #now process each observation (row in the table)
    for obs_vals,obs_id,obs_metadata in table.iterObservations():
        row=[]
        if max_len_metadata >0:
            row=process_metadata(obs_metadata[metadata_name],metadata_name,obs_id)
        
        #Add 'Unclassified' if the metadata doesn't fill each level
        len_defined_metadata=len(row)
        if len_defined_metadata < max_len_metadata:
            for i in range(max_len_metadata - len_defined_metadata):
                row.append('Unclassified')

        if include_obs_id:
            #Add the observation id as the last "Level"
            if obs_id.isdigit():
                #Need to add something to the id if it a number identfier (e.g. gg OTU ids)
                row.append('ID'+'_'+obs_id)
            else:
                row.append(obs_id)

        #Add count data to the row
        row.extend(map(str,obs_vals))
        print "\t".join(row)
        
    

if __name__ == "__main__":
    main()

