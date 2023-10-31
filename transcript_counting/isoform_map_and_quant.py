#!/usr/bin/env python3
'''Script for mapping individual cell fasta files to a ref transcriptome and quantifying isoforms with salmon'''
'''dependencies:minimap2/salmon'''



import argparse
import editdistance
import numpy as np
import os
import sys
import shutil
import mappy as mm
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to individual cell .final.fasta UMI-merged/polished files
parser.add_argument('-o', '--output_path', type=str) #where you want 'cellsams' directory to go 
parser.add_argument('-c', '--config_file', type=str) #config file containing paths to executables if not in $PATH
parser.add_argument('-r', '--ref_tx', type=str) #reference TRANSCRIPTOME for minimap2 alignment


args = parser.parse_args()
input_path=args.input_path
path=args.output_path
config_file=args.config_file
ref_tx=args.ref_tx

'''Parse config'''
def configReader(configIn):
    progs = {}
    for line in open(configIn): 
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t') 
        progs[line[0]] = line[1] 
        
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'samtools','featureCounts', 'psl2pslx']) #don't think I actually need psl2pslx and emtrey...
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
        #sys.stderr.write('Using ' + str(missing)
                         #+ ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
minimap2 = progs['minimap2']


'''Read each cell fasta in input path, create temp trimmed fasta for minimap2 alignment'''

def trim (inFile):
    readdict={}
    filename=inFile.split('/')
    filename=filename[len(filename)-1]
    cell=filename.split('_')[0]
    number=filename.split('_')[1]
    barcode=filename.split('_')[2]
    bc=barcode.split('.')[0]
    for line in open(inFile):
        if line.startswith('>'):
            name=line.split('_')[0]
            name=name+'_'+cell+'_'+number+'_'+bc
        else:
            sequence=line
            trim=len(sequence)-29
            sequence=sequence[0:trim]
            readdict[name]=sequence
    return readdict


def map_reads(name,inFile,ref_tx,path):
    #sys.stderr.write('Mapping reads to %s' %(ref_fasta))
    out= path + 'cellsams_tx/'+ name +'.sam'
    os.system('%s -ax map-ont -N 100 -p .99 %s %s > %s' % (minimap2, ref_tx, inFile, out))
    salmon_quant(out,ref_tx,path,name)

def salmon_quant(infile, ref_tx, path, name):
    os.system('salmon quant -t %s -l U -a %s --ont -o %s' % (ref_tx, infile, path))
    filename=path+'/quant.sf'
    shutil.move(filename,path+'/'+name+'.sf') 

def concat_fasta(input_path, path, ref_tx):
    #sys.stderr.write('Concatenating reads')
    fileList=[]
    for file in os.listdir(input_path):
        #final=open(path+'/temptrimmed.fasta','w')
        #if '.final.fasta' in file:
            #fileList.append(file)
    #print(sorted(fileList, key=lambda x: int(x.split('_')[1])))
    #for file in sorted(fileList, key=lambda x: int(x.split('_')[1]))
        name=file.split('.')[0]
        final=open(path+'/temptrimmed.fasta','w')
        temp=path+'/temptrimmed.fasta'
        fasta_file = input_path + '/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i] 
            final.write('%s\n%s\n' % (i,sequence))
        #final.close()
        map_reads(name,temp,ref_tx,path)
    #return path+'/allfinalreads.fasta'


def main(input_path, path, ref_tx):
    os.mkdir(path+'cellsams_tx/')

    catfile=concat_fasta(input_path, path, ref_tx)

main(input_path, path, ref_tx)
