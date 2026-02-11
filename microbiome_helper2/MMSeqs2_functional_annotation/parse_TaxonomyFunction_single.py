import numpy as np
import pandas as pd
import csv
from collections import defaultdict
#import _pickle as cPickle
#import pickle  as cPickle

try:
    import _pickle as cPickle
except ImportError:
    import pickle
    
import bz2
from operator import itemgetter

import argparse, sys, textwrap

parser=argparse.ArgumentParser()

parser.add_argument('--taxafile', help = 'File mapping the reads to Taxa')
parser.add_argument('--taxafiletype', help = 'Source program for Taxa identification, either one of kraken2 or megan')
parser.add_argument('--funcfile', help = 'File mapping reads to functional categories')
parser.add_argument('--funcfiletype', help = 'Source program for Function identification, either one of uniref, refseq, megan or COG')
parser.add_argument('--m8file', help = 'BLAST or BLAST-like tab delimited output file mapping reads to Refseq or uniref IDs; required for Normalization to RPKG')
parser.add_argument('--database', help = 'Path to the folder containing the database used for analysis, including any .pbz2 files with gene lengths etc')


parser.add_argument('--multisample', help = textwrap.dedent('''For running multiple samples at a time, input a text file; overrides the above 5 arguments
    Format of the text file should contain 6 columns
    
    '''))

parser.add_argument('--outputf', help = 'Start of the file name for writing all output files. Note that multiple will be saved with different file endings')

# parser.add_argument('--unstratified', help = 'Boolean Y|N to output unstratified metabolic functions; must choose at most one of --stratified or --unstratified')
# parser.add_argument('--stratified', help = 'Boolean Y|N to output metabolic functions stratified by taxa; must choose at most one of --stratified or --unstratified')
# parser.add_argument('--map2EC', help = 'Boolean Y|N to output a matrix of EC numbers in samples')



def main():
    args = parser.parse_args()
    
    multi = args.multisample
    # strat = args.stratified
    # unstrat = args.unstratified
    outfile = args.outputf
    database = args.database

    # map2ECflag = args.map2EC
    
    genlenf = bz2.BZ2File(database+'/GeneLength.pbz2', 'rb')
    genelendict = cPickle.load(genlenf)
    
    # RSgenlenf = bz2.BZ2File('RefGeneLength.pbz2', 'rb')
    # RSgenelendict = cPickle.load(RSgenlenf)
    # 
    # UPgenlenf = bz2.BZ2File('UnirefGeneLength.pbz2', 'rb')
    # UPgenelendict = cPickle.load(UPgenlenf)
    # 
    # COGgenlenf = bz2.BZ2File('COGDBGeneLength.pbz2', 'rb')
    # COGgenelendict = cPickle.load(COGgenlenf)
    
    ECmapf = bz2.BZ2File(database+'/ECmapped.pbz2', 'rb')
    Seq2ECdict = cPickle.load(ECmapf)
    
    # if (map2ECflag == "Y"):
    #     
    #     ECmapf = bz2.BZ2File(database+'/ECmapped.pbz2', 'rb')
    #     Seq2ECdict = cPickle.load(ECmapf)
    #     
    #     # RS2ECmapf = bz2.BZ2File('RsECmapped.pbz2', 'rb')
    #     # RS2ECdict = cPickle.load(RS2ECmapf)
    #     # 
    #     # UP2ECmapf = bz2.BZ2File('UpECmapped.pbz2', 'rb')
    #     # UP2ECdict = cPickle.load(UP2ECmapf)
    #     # 
    #     # COG2ECmapf = bz2.BZ2File('COGECmapped.pbz2', 'rb')
    #     # COG2ECdict = cPickle.load(COG2ECmapf)
    # else:
    #     map2ECflag = "N"
        
    
    # MCreport = args.MicrobeCensusReport
    # GEdict = {}
    # 
    # if (MCreport):
    #     GEdict = parseMicrobeCensusReport(MCreport)
    #     GEdictflag = "Y"
    #     print (GEdict)
    # else:
    #     print ("MicrobeCensus report not provided; RPKG will not be calculated; will normalize to RPKM")
    #     GEdictflag = "N"
    
    
    if not multi:
        taxadict = {}
        funcdict = {}

        taxafile = args.taxafile
        taxafiletype = args.taxafiletype
        funcfile = args.funcfile
        funcfiletype = args.funcfiletype
        m8file = args.m8file
        
        print ("Running single sample:", taxafile,taxafiletype,funcfile,funcfiletype)
        
        (taxadict,funcdict) = coreRun(taxafile,taxafiletype,funcfile,funcfiletype)
        
        print ("Total reads mapped to taxa: " + str(len(taxadict)))
        print ("Total reads mapped to functions: " + str(len(funcdict)))
        
        combinedDict = mergeTaxaFunc(taxadict,funcdict)
        print ("Reads mapped to either taxa OR functions: " + str(len(combinedDict)))
        filteredDict = dict(filter(lambda x: len(x[1]) == 2, combinedDict.items()))
                
        print ("Reads mapped to both taxa AND functions: " + str(len(filteredDict)))
 
        funchash = defaultdict(list)
        for key, value in sorted(filteredDict.items()):
            funchash[value[1]].append(key)
                
    else:
        print ("Running multiple samples from file:", multi)
        
        with open(multi, newline = '') as multif:                                                                                          
            multi_reader = csv.reader(multif, delimiter='\t')
            SampleFuncRpkgdict_strat, SampleFuncRpkgdict_unstrat, SampleFuncRpkgdict_strat_ec, SampleFuncRpkgdict_unstrat_ec = {}, {}, {}, {}
            for line in multi_reader:
                
                line_elements = len(line)
		
                print ("line:",line)
                print ("Total elements:",line_elements)
                if (line_elements != 6):
                    continue

                taxadict = {}
                funcdict = {}
                tot_reads_mapped,taxa_AND_func_mapped,taxa_mapped,taxa_OR_func_mapped,func_mapped,EC_mapped = 0,0,0,0,0,0
                perc_EC_mapped,perc_func_mapped,perc_taxa_AND_func_mapped,perc_taxa_mapped,perc_taxa_OR_func_mapped = 0,0,0,0,0
                sampletag = line[0]
                taxafile = line[1]
                taxafiletype = line[2]
                funcfile = line[3]
                funcfiletype = line[4]
                m8file = line[5]
                
                # if (funcfiletype == 'megan' and map2ECflag == "Y"):
                #     #genlenf = bz2.BZ2File('RefGeneLength.pbz2', 'rb')
                #     #genelendict = cPickle.load(genlenf)
                #     genelendict = RSgenelendict
                #     #RS2ECmapf = bz2.BZ2File('RsECmapped.pbz2', 'rb')
                #     #Seq2ECdict = cPickle.load(RS2ECmapf)
                #     Seq2ECdict = RS2ECdict
                # elif (funcfiletype == 'megan' and map2ECflag == "N"):
                #     #genlenf = bz2.BZ2File('RefGeneLength.pbz2', 'rb')
                #     #genelendict = cPickle.load(genlenf)
                #     genelendict = RSgenelendict
                # elif (funcfiletype == 'uniref' and map2ECflag == "Y"):
                #     #genlenf = bz2.BZ2File('UnirefGeneLength.pbz2', 'rb')
                #     #genelendict = cPickle.load(genlenf)
                #     genelendict = UPgenelendict
                #     #RS2ECmapf = bz2.BZ2File('UpECmapped.pbz2', 'rb')
                #     #RS2ECdict = cPickle.load(RS2ECmapf)
                #     Seq2ECdict = UP2ECdict
                # elif (funcfiletype == 'uniref' and map2ECflag == "N"):
                #     #genlenf = bz2.BZ2File('UnirefGeneLength.pbz2', 'rb')
                #     #genelendict = cPickle.load(genlenf)
                #     genelendict = UPgenelendict
                # elif (funcfiletype == 'COG' and map2ECflag == "Y"):
                #     genelendict = COGgenelendict
                #     Seq2ECdict = COG2ECdict
                # elif (funcfiletype == 'COG' and map2ECflag == "N"):
                #     genelendict = COGgenelendict
                # elif (funcfiletype == 'refseq' and map2ECflag == "Y"):
                #     genelendict = RSgenelendict
                #     Seq2ECdict = RS2ECdict
                # elif (funcfiletype == 'refseq' and map2ECflag == "N"):
                #     genelendict = RSgenelendict
                    
                
                print ("Current sample:", taxafile,taxafiletype,funcfile,funcfiletype,m8file)
                (taxadict,funcdict,genedict) = coreRun(taxafile,taxafiletype,funcfile,funcfiletype,m8file)
                
                tot_reads_mapped = len(genedict)
                
                print ("Total reads mapped: " + str(tot_reads_mapped))
                
                taxa_mapped = len(taxadict)
                func_mapped = len(funcdict)
                perc_taxa_mapped = (taxa_mapped/tot_reads_mapped)*100
                perc_func_mapped = (func_mapped/tot_reads_mapped)*100
                
                print ("Total reads mapped to taxa: " + str(len(taxadict)) + " percent " + str(perc_taxa_mapped))
                print ("Total reads mapped to functions: " + str(len(funcdict)) + " percent " + str(perc_func_mapped))
        
                combinedDict = mergeTaxaFunc(taxadict,funcdict)
                taxa_OR_func_mapped = len(combinedDict)
                perc_taxa_OR_func_mapped = (taxa_OR_func_mapped/tot_reads_mapped)*100
                
                print ("Reads mapped to either taxa OR functions: " + str(len(combinedDict)) + " percent " + str(perc_taxa_OR_func_mapped))
                
                filteredDict = dict(filter(lambda x: len(x[1]) == 2, combinedDict.items()))
                taxa_AND_func_mapped = len(filteredDict)
                perc_taxa_AND_func_mapped = (taxa_AND_func_mapped/tot_reads_mapped)*100

                print ("Reads mapped to both taxa AND functions: " + str(len(filteredDict)) + " percent " + str(perc_taxa_AND_func_mapped))
                
                first10pairs = {k: filteredDict[k] for k in list(filteredDict)[100:120]}
                print("resultant dictionary taxa and func: \n", first10pairs)

                
                funchash = defaultdict(list)
                for key, value in sorted(filteredDict.items()):
                    key = str(key).rstrip('\n')
                    val = str(value[1]).rstrip('\n')
                    funchash[val].append(key)
                            
                #first10pairs = {k: funchash[k] for k in list(funchash)[100:120]}
                #print("resultant dictionary : \n", first10pairs)
            
                print("Total unique functions: ", len(funchash))
                
                ge = len(filteredDict)
                ge = ge/1000000
                
                EChash = mapRS2EC(funchash,Seq2ECdict,genedict)
                    
                EC_mapped = len(EChash)
                perc_EC_mapped = (EC_mapped/tot_reads_mapped)*100
                    
                print ("Total ECs detected: ",str(len(EChash)) + " percent " + str(perc_EC_mapped))
                
                #SampleFuncRpkgdict_strat
                functaxahash = stratified(funchash,filteredDict)
                FuncRpkgDict = normalize_rpkg(functaxahash,genedict,genelendict,ge)
                SampleFuncRpkgdict_strat[sampletag] = FuncRpkgDict
                
                #SampleFuncRpkgdict_unstrat
                FuncRpkgDict = normalize_rpkg(funchash,genedict,genelendict,ge)
                SampleFuncRpkgdict_unstrat[sampletag] = FuncRpkgDict
                
                #SampleFuncRpkgdict_strat_ec
                ECtaxahash = stratified(EChash,filteredDict)
                FuncRpkgDict = normalize_rpkg(ECtaxahash,genedict,genelendict,ge)
                SampleFuncRpkgdict_strat_ec[sampletag] = FuncRpkgDict
                
                #SampleFuncRpkgdict_unstrat_ec
                FuncRpkgDict = normalize_rpkg(EChash,genedict,genelendict,ge)
                SampleFuncRpkgdict_unstrat_ec[sampletag] = FuncRpkgDict
                                
                # if (GEdictflag == "Y" and unstrat == "Y" and map2ECflag == "N"):
                #     ge = GEdict.get(sampletag)
                #     FuncRpkgDict = normalize_rpkg(funchash,genedict,genelendict,ge)
                #     SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # elif (GEdictflag == "Y" and strat == "Y" and map2ECflag == "N"):
                #       ge = GEdict.get(sampletag)
                #       functaxahash = stratified(funchash,filteredDict)
                #       FuncRpkgDict = normalize_rpkg(functaxahash,genedict,genelendict,ge)
                #       SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # elif (GEdictflag == "Y" and unstrat == "Y" and map2ECflag == "Y"):
                #     ge = GEdict.get(sampletag)
                #     EChash = mapRS2EC(funchash,Seq2ECdict,genedict)
                #     
                #     EC_mapped = len(EChash)
                #     perc_EC_mapped = (EC_mapped/tot_reads_mapped)*100
                #     
                #     print ("Total ECs detected: ",str(len(EChash)) + " percent " + str(perc_EC_mapped))
                #     #first10pairs = {k: EChash[k] for k in list(EChash)[10:20]}
                #     #print("resultant dictionary EC read hash: \n", first10pairs)
                #     FuncRpkgDict = normalize_rpkg(EChash,genedict,genelendict,ge)
                #     SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # elif (GEdictflag == "Y" and strat == "Y" and map2ECflag == "Y"):
                #     ge = GEdict.get(sampletag)
                #     EChash = mapRS2EC(funchash,Seq2ECdict,genedict)
                #     
                #     EC_mapped = len(EChash)
                #     perc_EC_mapped = (EC_mapped/tot_reads_mapped)*100
                #     
                #     print ("Total ECs detected: ",str(len(EChash)) + " percent " + str(perc_EC_mapped))
                #     
                #     ECtaxahash = stratified(EChash,filteredDict)
                #     FuncRpkgDict = normalize_rpkg(ECtaxahash,genedict,genelendict,ge)
                #     SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # elif (GEdictflag == "N" and unstrat == "Y" and map2ECflag == "N"):
                #     ge = len(filteredDict)
                #     ge = ge/1000000
                #     FuncRpkgDict = normalize_rpkg(funchash,genedict,genelendict,ge)
                #     SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # elif (GEdictflag == "N" and strat == "Y" and map2ECflag == "N"):
                #     ge = len(filteredDict)
                #     ge = ge/1000000
                #     functaxahash = stratified(funchash,filteredDict)
                #     FuncRpkgDict = normalize_rpkg(functaxahash,genedict,genelendict,ge)
                #     SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # elif (GEdictflag == "N" and unstrat == "Y" and map2ECflag == "Y"):
                #     ge = len(filteredDict)
                #     ge = ge/1000000
                #     EChash = mapRS2EC(funchash,Seq2ECdict,genedict)
                #     
                #     EC_mapped = len(EChash)
                #     perc_EC_mapped = (EC_mapped/tot_reads_mapped)*100
                #     
                #     print ("Total ECs detected: ",str(len(EChash)) + " percent " + str(perc_EC_mapped))
                #     #first10pairs = {k: EChash[k] for k in list(EChash)[10:20]}
                #     #print("resultant dictionary EC read hash: \n", first10pairs)
                #     FuncRpkgDict = normalize_rpkg(EChash,genedict,genelendict,ge)
                #     SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # elif (GEdictflag == "N" and strat == "Y" and map2ECflag == "Y"):
                #     ge = len(filteredDict)
                #     ge = ge/1000000
                #     EChash = mapRS2EC(funchash,Seq2ECdict,genedict)
                #     
                #     EC_mapped = len(EChash)
                #     perc_EC_mapped = (EC_mapped/tot_reads_mapped)*100
                #     
                #     print ("Total ECs detected: ",str(len(EChash)) + " percent " + str(perc_EC_mapped))
                # 
                #     ECtaxahash = stratified(EChash,filteredDict)
                #     FuncRpkgDict = normalize_rpkg(ECtaxahash,genedict,genelendict,ge)
                #     SampleFuncRpkgdict[sampletag] = FuncRpkgDict
                # #else:
                #     #print ("MicrobeCensus result not found; Will not normalize")
            
    
    all_dicts = [SampleFuncRpkgdict_strat, SampleFuncRpkgdict_unstrat, SampleFuncRpkgdict_strat_ec, SampleFuncRpkgdict_unstrat_ec]
    names = ['-strat-matrix-RPKM.txt', '-unstrat-matrix-RPKM.txt', '-strat-matrix-RPKM-withEC.txt', '-unstrat-matrix-RPKM-withEC.txt']
    
    print('Saving files now')
    
    for s in range(len(all_dicts)):
      SampleFuncRpkgdict = all_dicts[s]
      pdDF = pd.DataFrame.from_dict(SampleFuncRpkgdict, orient='index')
      pdDF.fillna(0, inplace = True)
      pdDFT = pdDF.T
      pdDFT.columns.name = 'function'
      print (pdDFT.columns)
      pd.DataFrame.to_csv(pdDFT, path_or_buf=outfile+names[s], sep='\t', na_rep='', header=True, index=True, index_label='function', mode='w', lineterminator='\n')
      
    
    # pdDF = pd.DataFrame.from_dict(SampleFuncRpkgdict, orient='index')
    # pdDF.fillna(0, inplace = True)
    # pdDFT = pdDF.T
    # pdDFT.columns.name = 'function'
    # print (pdDFT.columns)
    # 
    # #if (strat == "Y"):
    # #    pdDFT[['function','sequence']] = pdDFT['function'].str.split("|",expand=True)
    # 
    # pd.DataFrame.to_csv(pdDFT, path_or_buf=outfile, sep='\t', na_rep='', header=True, index=True, index_label='function', mode='w', line_terminator='\n', escapechar=None, decimal='.')
        

def mapRS2EC(FuncReadHash,rs2ecdict,read2rsdict):
    
    ec2readhash = defaultdict(list)
    read2echash = {}
    
    for func, readarray in FuncReadHash.items():
        for read in readarray:
            rsidarray = read2rsdict[read]
            #print ("RefSeq Ids: ", rsidarray)
            ECarray = getECsforRSIds(rsidarray, rs2ecdict)
            #print ("EC lis: ", ECarray)
            if (ECarray):
                #print ("EC lis: ", ECarray)
                MFreqEC = rankECsforRead(ECarray)
                read2echash[read]=MFreqEC
    
    for readid, ec in read2echash.items():
        ec = "EC:" + str(ec)
        ec2readhash[ec].append(readid)
    
    
    return ec2readhash
        

def getECsforRSIds(rsidsarr,rs2echash):
    # get a list of ECs mapped to RefSeq IDs for a given read
    fullECarray=[]
    default = ""
    for rsid in rsidsarr:
        ecarray = rs2echash.get(rsid,default)
        for ec in ecarray:
            fullECarray.append(ec)
    
    return fullECarray
    
def rankECsforRead(ecarray):
    # get the most frequent EC in the list for each read
    
    unique, counts = np.unique(ecarray, return_counts=True)
    #print ("Unique + counts: ", unique, counts)
    ecCounts = {k: v for k, v in zip(unique, counts)}
    
    sorted_ecCounts = sorted(ecCounts.items(), key = itemgetter(1),reverse = True)
    
    #print ("EC counts: ", sorted_ecCounts, type(sorted_ecCounts))
    
    mfreqECcountTup = sorted_ecCounts[0]
    
    mfreEC = mfreqECcountTup[0]
    
    return mfreEC
    
    
def stratified(FuncReadHash,filtFuncTaxadict):
    # Stratify the function hash, 
    # i.e break down each function into taxonomic classes and
    # return a new hash with function-taxa pair keys
    
    taxahash = defaultdict(list)
    for key, value in sorted(filtFuncTaxadict.items()):
        taxahash[key] = value[0]
    
    stratified_hash = defaultdict(list)
    for func,readarray in FuncReadHash.items():
        for read in readarray:
            taxon = taxahash[read]
            functaxakey = str(func).rstrip('\n') + "|" + str(taxon).rstrip('\n')
            
            stratified_hash[functaxakey].append(read)
                           
    return stratified_hash
    
    

def normalize_rpkg(dict_reads_mapped_to_func,dict_genes_mapped_to_read,refseq_gene_len_dict,GE):
    FuncRPKGdict = defaultdict(list)
    
    for func,readarray in dict_reads_mapped_to_func.items():
                    FuncLengthArray = []
                    for read in readarray:
                        genes = dict_genes_mapped_to_read[read]
                        gene_lengths = list(map(lambda x: refseq_gene_len_dict.get(x), genes))
                        list(np.float64(gene_lengths))
                        read_avg_gene_length = Average(gene_lengths)
                        #print (read,avg_gene_length)
                        #FuncAvgGeneLengthdict[read].append(avg_gene_length)
                        FuncLengthArray.append(read_avg_gene_length)

                    FuncAvgGeneLength = Average(FuncLengthArray)

                    numreads = len(readarray)
                    lengthkb = FuncAvgGeneLength/1000
                    
                    if (GE and lengthkb):
                        rpkg = numreads/lengthkb/GE
                    else:
                        rpkg = float('NaN')

                    FuncRPKGdict[func]=float(rpkg)
                    #print ("Function: " + str(func) + " RPKG: " + str(rpkg))
    
    return FuncRPKGdict

def parseMicrobeCensusReport(MCRfile):
    d = {}
    with open(MCRfile) as f:
        for line in f:
            fields = line.split("\t")
            S_tag = fields[0]

            if S_tag == "SampleTag":
                continue
            
            S_GE = float(fields[3])
            d[str(S_tag)] = S_GE    

    return d

def Average(lst):
    avg=0

    if(len(lst) > 0):
        avg = float(sum(np.float64(lst)) / len(lst))
            
    return avg

            
def mutate_dict(f,d):
    new_d = {}
    # apply a function to all elements of a dict
    # to generate a new dict
    for k, v in d.items():
        new_d[k] = f(v)
    
    return new_d

def coreRun(taxaf,taxaft,funcf,funcft,m8):
    
    taxadict = {}
    funcdict = {}

    print ("In Core ... >>>:",taxaft,funcft)
    
    if (taxaft == 'megan'):
        print ("Taxatype:Megan\n")
        taxadict=parseMeganTaxafile(taxaf)
    elif (taxaft == 'kraken2'):
        print ("Taxatype:Kraken2\n")
        taxadict=parseKraken2Taxafile(taxaf)
    else:
        print ("Taxa type not recognised\n")
        
        
    if (funcft == 'megan'):
        print ("Functype:Megan\n")
        funcdict = parseMeganFuncfile(funcf)
    elif (funcft == 'uniref'):
        print ("Functype:uniref\n")
        funcdict = parseUnirefFuncfile(funcf)
    elif (funcft == 'COG'):
        print ("Functype:COG\n")
        funcdict = parseUnirefFuncfile(funcf)
    elif (funcft == 'refseq'):
        print ("Functype:refseq\n")
        funcdict = parseUnirefFuncfile(funcf)
    else:
        print ("Func type not recognised\n")
    
    if (m8):
        print ("m8 file given...processing")
        genedict = parseBlastm8(m8)
    
    return taxadict,funcdict,genedict
    

def mergeTaxaFunc(dict1,dict2):
    # Combine the values with same keys 
    result = defaultdict(list)
    # 
    for d in (dict1,dict2):
        for k,val in d.items():
            result[k].append(val)

            
    #first10pairs = {k: result[k] for k in list(result)[10:30]}
    #print("resultant dictionary : \n", first10pairs) 
    return result

def parseBlastm8(filename):
    d = defaultdict(list)
    print ("In m8 parse ... ")
    
    with open(filename) as f:
        for line in f:
            #(key, val) = line.split("\t")
            fields = line.split("\t")
            key = fields[0]
            val = fields[1] 
            #d[str(key)] = str(val)
            d[key].append(val)
            
    #first10pairs = {k: d[k] for k in list(d)[10:20]}
    #print("resultant dictionary : \n", first10pairs) 
    return d

    
def parseUnirefFuncfile(filename):
    
    d = {}
    print ("In Uniref Parsed func ... ")
    with open(filename) as f:
        for line in f:
            #fields = line.split("\t")
            (key, val) = line.split("\t")
            d[str(key)] = str(val)
    return d


def parseMeganFuncfile(filename):
    d = {}
    print ("In Megan func ... ")
    with open(filename) as f:
        for line in f:
            #fields = line.split("\t")
            (key, val) = line.split("\t")
            d[str(key)] = str(val)
    return d

def parseMeganTaxafile(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            d[str(key)] = str(val)
    return d        

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

if __name__ == "__main__":
    main();


    
