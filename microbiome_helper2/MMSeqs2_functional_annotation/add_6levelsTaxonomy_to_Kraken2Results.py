import argparse, sys, textwrap
import re
import csv
from collections import defaultdict
import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

import urllib

#from taxadb.taxid import TaxID
from ete3 import NCBITaxa
from ete3 import PhyloTree
#from ete3 import NCBITaxa
ncbi = NCBITaxa()
try:
    response = urllib.request.urlopen('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/')
    ncbi.update_taxonomy_database()
except urllib.error.URLError as e:
    logger.error('URLError: %s', e)


parser=argparse.ArgumentParser()

parser.add_argument('--taxafilelist', help = 'List of files mapping the reads to Taxa')
parser.add_argument('--outfile', help = 'File to write the taxonomy count matrix')

def main():
    args = parser.parse_args()
    multifilelist = args.taxafilelist
    outfile = args.outfile
    
    with open(multifilelist, newline = '') as multif:                                                                                          
            multi_reader = csv.reader(multif, delimiter='\t')
            SampleFuncRpkgdict = {}
            
            for line in multi_reader:
                print ("line:",line)
                sampletag = line[0]
                filename = line[1]
                #taxid = TaxID(dbtype='sqlite', dbname='/home/dhwani/databases/taxadb.sqlite')
                #tftype = args.taxafiletype
                taxadict=parseKraken2Taxafile(filename)
                
                
                rank_search_string = "phylum,class,order,family,genus,species"
                prefix_dict = {
                    "phylum": "p_",
                    "class": "c_",
                    "order": "o_",
                    "family": "f_",
                    "genus": "g_",
                    "species": "s_"
                    }
                
                taxadictWlineage, taxaIdDict = get_full_lineage(taxadict, rank_search_string, prefix_dict)
                first10pairs = {k: taxaIdDict[k] for k in list(taxaIdDict)[1:20]}
                print("resultant dictionary : \n", first10pairs)
                
                SampleFuncRpkgdict[sampletag] = taxadictWlineage
                
                
    pdDF = pd.DataFrame.from_dict(SampleFuncRpkgdict, orient='index')
    pdDF.fillna(0, inplace = True)
    pdDFT = pdDF.T
    pdDFT.columns.name = 'lineage'
    print (pdDFT.columns)
    
    #if (strat == "Y"):
    #    pdDFT[['function','sequence']] = pdDFT['function'].str.split("|",expand=True)
    #new_d = pd.Series(d)
    
    pdDFT["ID"] = pd.Series(taxaIdDict)
    
    
    first_col = pdDFT.pop('ID')
    count_tracker = 0 
    unique_first_col = []

    for i in first_col:
        new_id = str(count_tracker)+"_"+str(i)
        unique_first_col.append(new_id)
        count_tracker+=1
        
    print ("First col : ", unique_first_col)
    
    pdDFT_mod=pdDFT.reset_index()
    pdDFT_mod.insert(0, 'newID', unique_first_col)
    
    pdDFT_mod.set_index('newID', inplace=True)
    pdDFT_mod.columns.name = 'newID'
    #pdDFT.index.name = 'newID'
    
    
    #print (pdDFT.columns)
    
    pd.DataFrame.to_csv(pdDFT_mod, path_or_buf=outfile, sep='\t', na_rep='', header=True, index=True, mode='w', line_terminator='\n', escapechar=None, decimal='.')

    

def parseKraken2Taxafile(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            fields = line.split("\t")
            classified = fields[0]
            key = fields[1]
            val = fields[2]    
            #print (fields[0])
            if (classified == 'C'):
                d[str(key)] = str(val)
    return d

def get_full_lineage(taxdict, rank_str, pref_dict):
    
    lineage_dict = defaultdict(list)
    lineage_id_dict = {}
    lineage = ""
    count="1"
    
    for k,val in taxdict.items():
        complete_taxa_string=""
        #s = 'asdf=5;iwantthis123jasd'
        #tax_id = re.search('(taxid (.*))', val)
        #print(tax_id.group(1))
        start = '(taxid '
        end = ')'
        #s = 'asdf=5;iwantthis123jasd'
        #print (val[val.find(start)+len(start):val.rfind(end)])
        id = val[val.find(start)+len(start):val.rfind(end)]
        #print("Tax ID is:",id)
        if (id == '0'):
            id = '1'
        
        try:
            lineage = ncbi.get_lineage(id)
            #import _pickle as cPickle
        except ValueError:
            #print("Oops!  That was no valid number.  Try again...")
            pass
            #import pickle
            

        if (lineage):
            ranks = ncbi.get_rank(lineage)
            #print ("ranks :", ranks)
            names = ncbi.translate_to_names(lineage)
            #print ("names :", names)
            
            id_name_dict = {lineage[i]: names[i] for i in range(len(lineage))}
            #print ("Id-name dict:", id_name_dict)
            rank_str_list = list(rank_str.split(","))
            newDict = {key: value for (key, value) in ranks.items() if value in rank_str_list }
            #print ("filetered: ", newDict)
            
            inv_FiltDict = {v: k for k, v in newDict.items()}
            
            all_taxa=[]
            
            for tax_rank in rank_str_list:
                tax_prefix = pref_dict[tax_rank]
                #print ("prefix ..: ", tax_prefix)
                
                if (tax_rank in inv_FiltDict):
                    rank_tax_id = inv_FiltDict[tax_rank]
                else:
                    rank_tax_id = 'NA'
                
                if (rank_tax_id in id_name_dict):
                    rank_tax_name = id_name_dict[rank_tax_id]
                else:
                    rank_tax_name = 'NA'
                
                rank_tax_name_w_pref = tax_prefix + "_" + rank_tax_name
                all_taxa.append(rank_tax_name_w_pref)
                
                #complete_taxa_string += complete_taxa_string + rank_tax_name
                #print ("rank:\t",tax_rank)
                #print ("name:\t",rank_tax_name)
            complete_taxa_string = ';'.join(all_taxa)

            #print ("Taxa string :", complete_taxa_string)
            
            lineage_dict[complete_taxa_string].append(count)
            lineage_id_dict[complete_taxa_string] = id
        else:
            print ("lineage not found, check ID:", id)
    
    #my_dictionary = dict(map(lambda kv: (kv[0], f(kv[1])), my_dictionary.iteritems()))

    lineage_count_dict = dict(map(lambda x: (x[0], len(x[1])), lineage_dict.items()))
    return lineage_count_dict, lineage_id_dict
    
if __name__ == "__main__":
    main();
