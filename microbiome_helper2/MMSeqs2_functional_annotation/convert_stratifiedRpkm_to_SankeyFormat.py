import argparse, sys, textwrap
import re
import csv
from collections import defaultdict
import numpy as np
import pandas as pd

#from taxadb.taxid import TaxID
from ete3 import NCBITaxa
from ete3 import PhyloTree
#from ete3 import NCBITaxa
ncbi = NCBITaxa()

parser=argparse.ArgumentParser()

parser.add_argument('--StratFileName', help = 'filename of file containing the stratified output')
#parser.add_argument('--filetype', help = 'type of files to parse, either of kraken2 OR stratified')
parser.add_argument('--taxaAbundFile', help = 'The krkaken2 taxa abundance file generated with 6 levels of taxonomy')
parser.add_argument('--outfile', help = 'File to write the taxonomy count matrix')

def main():
    args = parser.parse_args()
    filename = args.StratFileName
    #filetype = args.filetype
    taxafile = args.taxaAbundFile
    outfile = args.outfile
    
    #if (filetype == 'kraken2'):
    
        #with open(multifilelist, newline = '') as multif:                                                                                          
                #multi_reader = csv.reader(multif, delimiter='\t')
                #SampleFuncRpkgdict = {}
                
                #for line in multi_reader:
                    #print ("line:",line)
                    #sampletag = line[0]
                    #filename = line[1]
                    ##taxid = TaxID(dbtype='sqlite', dbname='/home/dhwani/databases/taxadb.sqlite')
                    ##tftype = args.taxafiletype
                    #taxadict=parseKraken2Taxafile(filename)
                    
                    
                    #rank_search_string = "phylum,class,order,family,genus,species"
                    #prefix_dict = {
                        #"phylum": "p_",
                        #"class": "c_",
                        #"order": "o_",
                        #"family": "f_",
                        #"genus": "g_",
                        #"species": "s_"
                        #}
                    
                    #taxadictWlineage, taxaIdDict = get_full_lineage(taxadict, rank_search_string, prefix_dict)
                    #first10pairs = {k: taxaIdDict[k] for k in list(taxaIdDict)[1:20]}
                    #print("resultant dictionary : \n", first10pairs)
                    
                    #SampleFuncRpkgdict[sampletag] = taxadictWlineage
                    
                    
        #pdDF = pd.DataFrame.from_dict(SampleFuncRpkgdict, orient='index')
        #pdDF.fillna(0, inplace = True)
        #pdDFT = pdDF.T
        #pdDFT.columns.name = 'lineage'
        #print (pdDFT.columns)
        
        ##if (strat == "Y"):
        ##    pdDFT[['function','sequence']] = pdDFT['function'].str.split("|",expand=True)
        ##new_d = pd.Series(d)
        
        #pdDFT["ID"] = pd.Series(taxaIdDict)
        
        
        #first_col = pdDFT.pop('ID')
        #count_tracker = 0 
        #unique_first_col = []

        #for i in first_col:
            #new_id = str(count_tracker)+"_"+str(i)
            #unique_first_col.append(new_id)
            #count_tracker+=1
            
        #print ("First col : ", unique_first_col)
        
        #pdDFT_mod=pdDFT.reset_index()
        #pdDFT_mod.insert(0, 'newID', unique_first_col)
        
        #pdDFT_mod.set_index('newID', inplace=True)
        #pdDFT_mod.columns.name = 'newID'
        ##pdDFT.index.name = 'newID'
        
        
        ##print (pdDFT.columns)
        
        #pd.DataFrame.to_csv(pdDFT_mod, path_or_buf=outfile, sep='\t', na_rep='', header=True, index=True, mode='w', line_terminator='\n', escapechar=None, decimal='.')

    #elif (filetype == 'stratified'):
    keyStrDict = {}
    idFulltaxaDict = {}
    lineage_otuidDict = {}
    genecopynoDict = {}
    
    #filename = multifilelist
    
    DF = pd.read_table(filename, index_col=0)
    
    taxaAbundDF = parseTaxaAbundancefile(taxafile)
    
    #idFulltaxaDict = pd.Series(taxaAbundDF.Letter.values,index=df.Position).to_dict()
    idFulltaxaDict = dict(zip(taxaAbundDF['index'], list(taxaAbundDF.index)))
    
    first10pairs = {k: idFulltaxaDict[k] for k in list(idFulltaxaDict)[1:20]}
    print("otu id - lineage dict dictionary : \n", first10pairs)
    #print ("element of taxaDF: ", taxaAbundDF.at["6_6072","BB16-15D-2_mmseqs"])
    #print ("indexex of taxaDf: ", list(taxaAbundDF.index.values))
    
    

    f = open(filename, "r")
    first_line = f.readline().strip('\n')
    for line in f:
        fields = line.split("\t")
        FuncTaxastr = fields[0]
        FuncTaxa = FuncTaxastr.split("|")
        func = FuncTaxa[0]
        TaxaStr = FuncTaxa[1]
        
        start = '(taxid '
        end = ')'
        #s = 'asdf=5;iwantthis123jasd'
        #print (val[val.find(start)+len(start):val.rfind(end)])
        id = TaxaStr[TaxaStr.find(start)+len(start):TaxaStr.rfind(end)]
        if (id == '0'):
            id = '1'
        
        rank_search_string = "phylum,class,order,family,genus,species"
        prefix_dict = {
                    "phylum": "p_",
                    "class": "c_",
                    "order": "o_",
                    "family": "f_",
                    "genus": "g_",
                    "species": "s_"
                    }
        #print ("Function\tTaxa\tTaxid", func,TaxaStr,id)
        replacementStr_compLineage = get_full_lineage_for_id(id,rank_search_string,prefix_dict)
        #print ("Replacement string :", replacementStr_compLineage)
        newFuncTaxaStr = func+"|"+replacementStr_compLineage
        #print ("new fucn taxa : ", newFuncTaxaStr)
        keyStrDict[FuncTaxastr] = newFuncTaxaStr
    
    first10pairs = {k: keyStrDict[k] for k in list(keyStrDict)[1:20]}
    print("orig id  - lineage dictionary : \n", first10pairs)
    #lineageOtuId_dict = {value: key for key, value in .items()}
    #Insert the new key string into df and remove the old key string
    
    for orig_id,new_id in keyStrDict.items():
        func_taxa = new_id.split("|")
        lineage_str = func_taxa[1]
        #otu_id = idFulltaxaDict[lineage_str]
        otu_id = idFulltaxaDict.get(lineage_str)
        if otu_id:
            lineage_otuidDict[orig_id] = otu_id
            genecopynoDict[orig_id] = 1.0
        else:
            print ("Damn OTU id not found")
        
    
    first10pairs = {k: lineage_otuidDict[k] for k in list(lineage_otuidDict)[1:20]}
    print("orig id  - lineage dictionary : \n", first10pairs)


            #d2[key] = d1[key]

    DF["newID"] = pd.Series(keyStrDict)
    #DF["OTU"] = pd.Series(lineage_otuidDict)
    #DF["GeneCountPerGenome"] = pd.Series(genecopynoDict)
    
    col_name="newID"
    first_col = DF.pop(col_name)
    #DF=DF.reset_index()
    
    #print ("new id column :\n", first_col)
    
    
    DF.insert(0, col_name, first_col)
    
    #DF_mod.set_index(col_name, inplace=True)
    #DF_mod.columns.name = col_name
    DF=DF.reset_index()
    DF = DF.drop('function', 1)
    
    #split first column into "function" and "taxa"
    
    DF[['Gene','Genus']] = DF['newID'].str.split('|',expand=True)
    DF = DF.drop('newID', 1)
    
    # rearrange format by melting sample and taxa
    #stubnames = sorted(set([match[0] for match in DF.columns.str.findall(r'[B]\(.*\)').values if match != []]))
    #l = pd.wide_to_long(DF, stubnames='contributions', i=['famid', 'birth'], j='age')
    #print ("Stubnames: ", stubnames)
    
    DF_long = DF.melt(id_vars=["Gene", "Genus"], var_name = "Sample", value_name = "Contribution")
#        DF_long2 = DF.melt(id_vars=["Gene", "Genus"], var_name = "OTU", value_name = "OTU")
    
    
    cols = DF_long.columns.tolist()
    cols[0],cols[2] = cols[2],cols[0]
    print("Final Columns of DF:", cols)
    
    DF_long = DF_long[cols]
    
    #DF_long[~DF_long.Sample.str.contains("OTU", na=False)]
    
    pd.DataFrame.to_csv(DF_long, path_or_buf=outfile, sep='\t', na_rep='', header=True, index=False, mode='w', line_terminator='\n', escapechar=None, decimal='.')

def parseTaxaAbundancefile(taxaabundfile):
    taxDF = {}
    taxDF = pd.read_table(taxaabundfile, index_col = 'newID', header = 0)
    cols = taxDF.columns.tolist()
    print ("tax DF columns :", cols)
    #taxDF_long = taxDF.melt(id_vars=['newID', "index"], var_name="Sample", value_name="Abund")
    pd.DataFrame.to_csv(taxDF, path_or_buf='taxtest.out', sep='\t', na_rep='', header=True, index=True, mode='w', line_terminator='\n', escapechar=None, decimal='.')
    #taxDF = taxDF.reset_index
    #taxDF.set_index([taxDF.iloc[0], taxDF.columns[0]])
    #taxDF.set_index(list(taxDF)[0])
    #taxDF.set_index(taxDF.iloc[0].values)
    return taxDF
    #pd.DataFrame.to_csv(taxDF_long, path_or_buf='taxtest.out', sep='\t', na_rep='', header=True, index=False, mode='w', line_terminator='\n', escapechar=None, decimal='.')
    #return taxDF_long



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

def get_full_lineage_for_id(idStr, rank_str, pref_dict):
    lineage = ""
    try:
        lineage = ncbi.get_lineage(idStr)
            #import _pickle as cPickle
    except ValueError:
        print("Oops!  That was no valid number.  Try again...", idStr)
        pass
    
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
    return complete_taxa_string
    


def get_full_lineage(taxdict, rank_str, pref_dict):
    lineage = ""
    lineage_dict = defaultdict(list)
    lineage_id_dict = {}
    
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
        if (id == '0'):
            id = '1'
        
        try:
            lineage = ncbi.get_lineage(id)
            #import _pickle as cPickle
        except ValueError:
            print("Oops!  That was no valid number.  Try again...")
            pass
            #import pickle
            
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
    
    #my_dictionary = dict(map(lambda kv: (kv[0], f(kv[1])), my_dictionary.iteritems()))

    lineage_count_dict = dict(map(lambda x: (x[0], len(x[1])), lineage_dict.items()))
    return lineage_count_dict, lineage_id_dict
    
if __name__ == "__main__":
    main();
