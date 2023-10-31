#!/usr/bin/env python3
'''Take cell.fasta files from input dir, trim ends, cat, and map to transcriptome'''


import argparse
#import editdistance
#import numpy as np
import os
import sys
#import mappy as mm
#import shutil
#from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to cell-type directory
parser.add_argument('-o', '--output_path', type=str) #where you want temp output files and assignedreads directory to go
parser.add_argument('-r', '--ref_transcriptome', type=str) #reference GENOME for splice-aware minimap2 alignment
parser.add_argument('-n', '--sample_name', type=str) #sample name

args = parser.parse_args()
#input_file=args.input_file
input_path=args.input_path
path=args.output_path
#config_file=args.config_file
trx=args.ref_transcriptome
name=args.sample_name


'''Parse each cell fasta, trim barcode and UMI sequences from each read and concat into one fasta file'''
def trim (inFile):
    readdict={}
    filename=inFile.split('/')
    filename=filename[len(filename)-1]
    #rank=filename.split('_')[0]
    #cell=filename.split('_')[1]
    #number=filename.split('_')[1]
    #barcode=filename.split('_')[2]
    #bc=barcode.split('.')[0]
    for line in open(inFile):
        if line.startswith('>'):
            name=line.split('_')[0]
        else:
            sequence=line
            trim=len(sequence)-29
            sequence=sequence[0:trim]
            readdict[name]=sequence
    return readdict

def concat_fasta(input_path, path, name):
    #sys.stderr.write('Concatenating reads')
    final = open(path + '/'+name+'.fasta', 'w')
    fileList=[]
    for file in os.listdir(input_path):
        if 'final.fasta' in file:
            fileList.append(file)
    for file in sorted(fileList, key=lambda x: (x)):
        fasta_file = input_path + '/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i]
            final.write('%s\n%s\n' % (i,sequence))
    return path+'/'+name+'.fasta'


'''Map reads in the concatenated fasta file with minimap2'''
def map_reads (inFile,trx,path,name):
    out= path + '/'+name+'.sam'
    os.system('%s -ax map-ont -N 100 -p .99\
              %s %s > %s 2> ./minimap2_messages.txt'
              % ('minimap2',trx, inFile, out))
    return out


'''Convert sam alignment file to bam'''
def sam_con (inFile):
    #sys.stderr.write('Converting sam to bam file')
    outname=inFile.split('.')[0]
    os.system('samtools view -S -b %s > %s.bam' %(inFile, outname))
    os.system('samtools sort -o %s.sorted.bam %s.bam' %(outname, outname))
    os.system('samtools index -b %s.sorted.bam' %(outname))

def main(input_path, path, trx,name):
    catfile=concat_fasta(input_path, path,name)
    #sys.stderr.write('Trimming and concatenating fastas')
    samfile=map_reads(catfile, trx, path,name)
    sam_con(samfile)
    #bamfile=path+'/'+'UC2_Neurons.bam'


main(input_path, path, trx, name)

