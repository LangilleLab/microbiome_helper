#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 21:41:11 2020

@author: Dhwani
"""
# import numpy as np
# import pandas as pd
import csv
from collections import defaultdict
# import _pickle as cPickle
# import bz2
from operator import itemgetter
import re

import argparse, sys, textwrap

parser=argparse.ArgumentParser()

parser.add_argument('--rpkfile', help = 'Stratified or Unstratified functional abundance File')
parser.add_argument('--protIdfile', help = 'COG file mapping acc to GI')
parser.add_argument('--Cogfile', help = 'COG DB csv File')
parser.add_argument('--CogNamesfile', help = 'COG file mapping COD Ids to names')
parser.add_argument('--outfilename', help = 'Output file name with added descriptions')


def main():
    args = parser.parse_args()
    
    rpkfilename = args.rpkfile
    outp = args.outfilename
    
    
    protIdfilename = args.protIdfile
    
    cogfilename = args.Cogfile
    
    cognamesfilename = args.CogNamesfile
    
    cogProtIdhash = parseProtIds(protIdfilename)
    
    #first10pairs = {k: cogProtIdhash[k] for k in list(cogProtIdhash)[10:20]}
    #print("resultant dictionary : \n", first10pairs) 
    
    cognameshash = parseCognames(cognamesfilename)
    
    cogDBhash = parseCOGDB(cogfilename)
    
    ofh = open(outp, "w+")
    default = ""
    with open(rpkfilename , 'r') as rpkf:
        first_line = rpkf.readline()
        ofh.write(first_line)
        for line in rpkf:
            fields = re.split(r'\t+', line.rstrip('\n'))
            full_id = fields[0]
            #print(id)
            m = bool(re.search("\|", full_id))
            # print (str(m))
            if m:
                # print("found stratified file ...")
                parts = re.split("\|", full_id)
                id = parts[0]
                tax = parts[1]
                id = re.sub('\..*', '', id)
                # print (len(parts),tax)
                protID = cogProtIdhash.get(id, default)
                # print (protID)
                cogID = cogDBhash.get(protID,default)
                cogname = cognameshash.get(str(cogID),default)
                # print (cogname)
                lst_to_print = (protID,cogID,cogname)
                newID = ('-'.join(lst_to_print))
                newID = newID + "|" + str(tax)
                # print (newID)
    
            else:
                # print ("found unstr")
                id = full_id
                id = re.sub('\..*', '', id)
                # print (len(parts),tax)
                protID = cogProtIdhash.get(id, default)
                # print (protID)
                cogID = cogDBhash.get(protID,default)
                cogname = cognameshash.get(str(cogID),default)
                # print (cogname)
                lst_to_print = (protID,cogID,cogname)
                newID = ('-'.join(lst_to_print))
                # newID = newID + "|" + str(tax)
                # print (newID)
            
            fields[0] = newID
            
            lst_to_print_final = ('\t'.join(fields))
            # print (lst_to_print_final)
            ofh.write(lst_to_print_final + '\n')
def parseCOGDB(filename):
    d = defaultdict(list)
    print ("In COG DB parse ... ")
    
    with open(filename) as f:
        for line in f:
            #(key, val) = line.split("\t")
            fields = line.split(",")
            key = str(fields[2]).rstrip('\n')
            val = str(fields[6]).rstrip('\n') 
            # key = fields[2]
            # val = fields[6]
            d[key] = val
            #d[str(key).rstrip('\n')]=str(val).rstrip('\n')
            
    first10pairs = {k: d[k] for k in list(d)[10:20]}
    print("resultant dictionary : \n", first10pairs) 
    return d


def parseProtIds(filename):
    d = defaultdict(list)
    print ("In COG protein parse ... ")
    
    with open(filename) as f:
        for line in f:
            #(key, val) = line.split("\t")
            fields = line.split("\t")
            key = str(fields[1]).rstrip('\n')
            val = str(fields[0]).rstrip('\n') 
            # key = fields[1]
            # val = fields[0]
            d[str(key)] = val
            # d[str(key).rstrip('\n')]=str(val).rstrip('\n')
            
    first10pairs = {k: d[k] for k in list(d)[10:20]}
    print("resultant dictionary protIDs: \n", first10pairs) 
    return d

    
def parseCognames(filen):
    d = defaultdict(list)
    print ("In COG names parse ... ")
    
    with open(filen) as fn:
        for line in fn:
            #(key, val) = line.split("\t")
            if not re.match("#",line):
                fields = line.split("\t")
                key = str(fields[0]).rstrip('\n')
                val = str(fields[2]).rstrip('\n') 
                # key = fields[0]
                # val = fields[2]
                d[key] = val
                # d[str(key).rstrip('\n')]=str(val).rstrip('\n')
            
    first10pairs = {k: d[k] for k in list(d)[10:20]}
    print("resultant dictionary : \n", first10pairs) 
    return d
    
    


if __name__ == "__main__":
    main();


