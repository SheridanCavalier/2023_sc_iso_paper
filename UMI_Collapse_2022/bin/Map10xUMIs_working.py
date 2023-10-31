#!/usr/bin/env python3

import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm
import shutil
import multiprocessing as mp 
from tqdm import tqdm

'''Evaluate UMIs for reads from the same cell assigned to the same gene feature if UMIs are L=<1, collapse and polish using the 'best' read as the template'''
def eval_UMIs(sorted_file, demuxed_file,path,racon):
    #basically I want to make a tuple of the gene dict where in each tuple the first value is the gene and the second value is a list of the UMIs
    cell=sorted_file.split('_')[1]
    barcode=(sorted_file.split('_')[2]).split('.')[0]
    finalfasta = open(path + '/'+'cell_'+ cell + '_'+barcode+'.final.fasta', 'a') #if you already have this file it will be APPENDED so if you are doing this more than once be mindful
    genedict={}
    previous_gene=''
    read_dict={}
    #maybe for the mp part we just give the gene name and the UMI and then the output is the UMI groups
    for line in open(sorted_file):
        line=line.rstrip()
        gene=line.split('\t')[1]
        name=line.split('\t')[0]
        sequence=line.split('\t')[2]
        read_dict[name]=sequence
        UMI=sequence[len(sequence)-17::-1][0:12] #for dropseq; UMI=sequence[len(sequence)-13::-1][0:8]
        if gene != previous_gene:
            previous_gene=gene
            genedict[gene]=[]
            genedict[gene].append(name)
            genedict[gene].append(sequence)
            genedict[gene].append(UMI)
        else:
            genedict[gene].append(name)
            genedict[gene].append(sequence)
            genedict[gene].append(UMI) #creates genedict that has gene_id keys storing the names, seqs, and UMIs of associated reads (in that order)
    new_dict={}
    for key in genedict: #key is referring to gene id
        length=len(genedict[key])
        if length < 6: #if there is only one read for this gene in the cell...
            finalfasta.write('>%s\n%s\n' % (genedict[key][0], genedict[key][1]))
        elif length >= 6:
            new_dict[key]=genedict[key]
    
    for entry in new_dict:
        make_tup(new_dict[entry],read_dict,path,sorted_file,demuxed_file,cell,barcode,racon)
    finalfasta.close()










def make_tup(entry,read_dict,path,sorted_file,demuxed_file,cell,barcode,racon): #need to fix the hashable error with this read_dict
    read_dict=read_dict
    finalfasta = open(path + '/'+'cell_'+ cell + '_'+barcode+'.final.fasta', 'a')
    group_dict={}
    tup=[]
    for i in range(0,len(entry),3):
        p=[(entry[i],entry[i+2])]
        tup=tup+p
    group_dict=sort_tuple_vals(tup)
    for UMI in group_dict:
        if len(group_dict[UMI])==1:
            readname,readseq=group_dict[UMI][0],read_dict[group_dict[UMI][0]]
            finalfasta.write('>%s\n%s\n' % (readname, readseq))

def sort_tuple_vals(my_tup):
    my_tup.sort(key = lambda x: x[1])
    dict={}
    whitelist=[]
    unmatched=[]
    previous=''
    for entry in my_tup:
        if previous=='':
            previous=entry[1]
            dict[entry[1]]=[entry[0],]
        elif entry[1]==previous:
            dict[entry[1]].append(entry[0])
        elif entry[1] != previous:
            dict[entry[1]]=[entry[0],]
            previous=entry[1]
    for key in dict:
        if len(dict[key])>1:
            whitelist.append(key)
        else:
            unmatched.append(key)
    #print(whitelist,unmatched) #whitelist,unmatched,dict)

    for i in range(0,len(whitelist)):
        for j in range(0,len(unmatched)):
            if editdistance.eval(unmatched[j],whitelist[i])== 1:
                #dict[whitelist[i]].append(dict[unmatched[j]])
                if unmatched[j] in dict.keys():
                    dict[whitelist[i]].append(dict[unmatched[j]])
                    del dict[unmatched[j]]
                else:
                    continue
    return dict


