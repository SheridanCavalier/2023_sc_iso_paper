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

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

from Map10xUMIs_working import eval_UMIs
#from make_tup import make_tup 


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--assigned_path', type=str) #path to /assignedreads folder that contains mapped and sorted fastas
parser.add_argument('-o', '--output_path', type=str) #where temp files and final.fasta files will go
parser.add_argument('-d', '--demuxed_path', type=str) #path to /demuxed/fastq files
parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH


args = parser.parse_args()
demuxed=args.demuxed_path
assigned=args.assigned_path
path=args.output_path
config_file=args.config_file

'''Parse config'''
def configReader(configIn):
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]

    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'samtools','featureCounts', 'psl2pslx']) #don't think I actually need psl2pslx and emtr$
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
        if key not in possible:
            raise Exception('Check config file')

    for missing in possible-inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
    return progs

progs = configReader(config_file)
minimap2 = progs['minimap2']
racon = progs['racon']
samtools = progs['samtools']
featureCounts = progs['featureCounts']


def iterate(assigned, demuxed, path):
    fileList1=[]
    fileList2=[]
    for file in os.listdir(assigned):
        if '.sorted' in file:
            fileList1.append(assigned+'/'+file)
    for file in os.listdir(demuxed):
        if '.fastq' in file:
            fileList2.append(demuxed+'/'+file) 
    fileList1= sorted(fileList1, key=lambda x: int(x.split('_')[1]))
    fileList2= sorted(fileList2, key=lambda x: int(x.split('_')[1]))

    for i in range(0,len(fileList2)):
        eval_UMIs(fileList1[i], fileList2[i],path, racon)
    #return mp_list,genedict

#def main(assigned,demuxed,path): #add mp module to the UMI collapse script because it's like the wrapper for everything
iterate(assigned,demuxed,path)
    #pool = mp.Pool(mp.cpu_count())
    #from process_genes import process_genes
    #with Pool() as pool:
        #pool.starmap(calculate, itertools.product(range(3), range(5), range(4), range(10)))
    #pool.starmap(process_genes,args = (mp_list,genedict))
    #result = pool.map(process_genes, list)


#if __name__ == '__main__':
    #main(assigned,demuxed,path) #get these values from the eval UMIs fxn and then make gene_processing it's own function and mp that
    #fxn(item,dict)
