from collections import defaultdict
from operator import itemgetter

import argparse, sys, textwrap
import os

parser=argparse.ArgumentParser()

parser.add_argument('--unirefm8Dir', help = 'Folder containing m8 files mapping the reads to Uniref')
parser.add_argument('--output_path', help = 'Path to the directory where the parsed output will be stored')

def main():
    args = parser.parse_args()
    #unirefFileList = args.unirefm8
    outpath = args.output_path
    folderN=args.unirefm8Dir
    #print (unirefFileList, type(unirefFileList))
    
    #for m8 in unirefFileList:
    for root, dirs, files in os.walk(folderN, topdown=False):
        
        for m8 in files:
            filepath=str(os.path.join(root, m8))
            outfilename = outpath + "/" + m8 + "-parsed.txt"
            parseunirefm8(filepath,outfilename)
        #first10pairs = {k: reads_to_UPids_hash[k] for k in list(reads_to_UPids_hash)[10:20]}
        #print("resultant dictionary readsToUP: \n", first10pairs)
        #for re in first10pairs:
            #tup_list = first10pairs[re]
            #print (type(tup_list))
    

def parseunirefm8(filename,outf):
    d = defaultdict(list)
    ofh= open(outf,"w+")
    #read_UP_dict = {}
    #print ("In m8 parse ... ")
    
    with open(filename) as f:
        for line in f:
            #(key, val) = line.split("\t")
            fields = line.split("\t")
            read = fields[0]
            UPid = fields[1]
            Evalue = fields[10]
            #d[str(key)] = str(val)
            d[read].append((UPid,Evalue))
            
    #first10pairs = {k: d[k] for k in list(d)[10:20]}
    #print("resultant dictionary : \n", first10pairs) 
    
    for re in d:
        tup_list = d[re]
        top_UP_id = pick_top_hit(tup_list)
        #first10pairs = {k: read_UP_Eval_hash[k] for k in list(read_UP_Eval_hash)[1:10]}
        #print(re, " resultant dictionary for read: \n", first10pairs)
        lst_to_print=(re,top_UP_id)
        ofh.write('\t'.join(lst_to_print) + '\n')
        #print (re,"\tTop hit: ", top_UP_id)
    
    return d

def pick_top_hit(tupList):
    di = {}
    di = dict(tupList)
    #print(type(di))
    mod_di = mutate_dict(lambda x: float(x), di)
    sorted_by_evalue = sorted(mod_di.items(), key = itemgetter(1))
    top_hit_tuple = sorted_by_evalue[0]
    top_hit_UP = top_hit_tuple[0]
    #first10pairs = {k: sorted_by_evalue[k] for k in list(sorted_by_evalue)[1:10]}
    #print(" resultant dictionary for read: \n", sorted_by_evalue)
    return top_hit_UP

def mutate_dict(f,d):
    new_d = {}
    # apply a function to all elements of a dict
    # to generate a new dict
    for k, v in d.items():
        new_d[k] = f(v)
    
    return new_d


if __name__ == "__main__":
    main();
