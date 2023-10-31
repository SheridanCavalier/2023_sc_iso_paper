#!/usr/bin/env python3

'''Take cell.fasta files from input dir, trim ends, cat, and map to genome'''

import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to cell-type directory
parser.add_argument('-o', '--output_path', type=str) #where you want temp output files and assignedreads directory to go
parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH
parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for splice-aware minimap2 alignment


args = parser.parse_args()
input_path=args.input_path
path=args.output_path
config_file=args.config_file
ref_genome=args.ref_genome

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
        #sys.stderr.write('Using ' + str(missing)
                         #+ ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
#poa = progs['poa']
minimap2 = progs['minimap2']
racon = progs['racon']
#consensus = progs['consensus']
samtools = progs['samtools']
featureCounts = progs['featureCounts']

'''Parse each cell fasta, trim barcode and UMI sequences from each read and concat into one fasta file'''
def trim (inFile):
    readdict={}
    filename=inFile.split('/')
    filename=filename[len(filename)-1]
    #print(filename)
    rank=filename.split('_')[0]
    cell=filename.split('_')[1]
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

def concat_fasta(input_path, path):
    #sys.stderr.write('Concatenating reads')
    final = open(path + '/naive.fasta', 'w')
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
    return path+'/naivefasta'


'''Map reads in the concatenated fasta file with minimap2'''
def map_reads (inFile,ref_fasta,path):
    #sys.stderr.write('Mapping reads to %s' %(ref_fasta))
    out= path + '/naive.sam'
    os.system('%s --secondary=no -ax splice \
              %s %s > %s 2> ./minimap2_messages.txt'
              % (minimap2,ref_fasta, inFile, out))
    return out


'''Convert sam alignment file to bam'''
def sam_con (inFile):
    #sys.stderr.write('Converting sam to bam file')
    outname=inFile.split('.')[0]
    os.system('samtools view -S -b %s > %s.bam' %(inFile, outname))
    os.system('samtools sort -o %s.sorted.bam %s.bam' %(outname, outname))
    os.system('samtools index -b %s.sorted.bam' %(outname))

def main(input_path, path, ref_genome):
    catfile=concat_fasta(input_path, path)
    #sys.stderr.write('Trimming and concatenating fastas')
    samfile=map_reads(catfile, ref_genome, path)
    sam_con(samfile)
    bamfile=path+'/'+'naive.bam'


main(input_path, path, ref_genome)

