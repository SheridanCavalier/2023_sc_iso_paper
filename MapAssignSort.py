#!/usr/bin/env python3

import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm
import shutil 
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to individual cell .fasta and _subs.fastq output files from Vollmers demultiplex_kmer script
                                                    #final.fasta files with mapped and collapsed UMIs will go into this folder as well
parser.add_argument('-o', '--output_path', type=str) #where you want temp output files and assignedreads directory to go 
parser.add_argument('-g', '--gtf_file', type=str) #path to reference genome gtf file for FeatureCount gene id assignments
parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH
parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for splice-aware minimap2 alignment


args = parser.parse_args()
#input_file=args.input_file
input_path=args.input_path
path=args.output_path
config_file=args.config_file
ref_genome=args.ref_genome
gtf=args.gtf_file


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
#poa = progs['poa'] 
minimap2 = progs['minimap2']
#racon = progs['racon']
#consensus = progs['consensus']
samtools = progs['samtools']
featureCounts = progs['featureCounts']


'''Parse each cell fasta, trim barcode and UMI sequences from each read and concat into one fasta file'''
def trim (inFile):
    readdict={}
    filename=inFile.split('/')
    filename=filename[len(filename)-1]
    #print(filename)
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


def concat_fasta(input_path, path):
    #sys.stderr.write('Concatenating reads')
    final = open(path + '/allreads.fasta', 'w')
    fileList=[]
    for file in os.listdir(input_path):
        if '.fasta' in file and 'final' not in file:
            fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):
        fasta_file = input_path + '/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i] 
            final.write('%s\n%s\n' % (i,sequence))
    return path+'/allreads.fasta'


'''Map reads in the concatenated fasta file with minimap2'''
def map_reads (inFile,ref_fasta,path):
    #sys.stderr.write('Mapping reads to %s' %(ref_fasta))
    out= path + '/allreads.sam'
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


'''Assign mapped reads to gene names with featureCount'''
def fcount (inFile, gtf):
    #sys.stderr.write('Getting gene assignments for mapped reads')
    os.system('featureCounts --primary -R CORE --maxMOp 1000 -Q 0 -L -t exon -g gene_id -a %s -o mysample_featureCount.txt %s' %(gtf, inFile))
#-L denotes long-reads, -Q 9 only assigns alignments with a MAPQ score greater than 9

'''Group assigned reads by cell barcode and gene assignment'''
def group_reads (inFile, input_path, output_path):
    #if not os. path. isdir(path+'assignedreads'):
        #shutil.rmtree(path+'assignedreads')
    #else:
    os.mkdir(path+'assignedreads')
 
    previous_cell=''
    dict={}
    ass_dict={} #haha
    fileList=[]
    bcdict={}
    fullnamedict={}
    for line in open(inFile): #inFile is allreads.bam.featurecounts
        name=line.rstrip().split('\t')[0] #full name
        rname=name.split('_')[0]
        if line.rstrip().split('\t')[3]!='NA'and rname not in ass_dict: #just take the primary assignment
            ass_dict[rname]=line.rstrip().split('\t')[3] #assignment dictionary at that fullname stores the gene id for that read--reads are being duplicated! Need to figure out how to fix this!
        cell=name.split('_')[2] #takes cell number from each read name (doesn't care about assignment)
        bc=name.split('_')[3] #takes the barcode from each read name (doesn't care about assignment)
        bcdict[cell]=bc #stores barcode information under cell number in the barcode dictionary
        #final = open(path + '/' + cell + '.new.fasta', 'w')
        if int(line.rstrip().split('\t')[2])==1:
            if cell!=previous_cell:
                previous_cell=cell
                dict[cell]=[]
                bcdict[cell]=bc
                if name not in dict.get(cell): #my clumsy way of preventing read duplicates and taking the first gene_id in featurecounts (could be better)
                    dict[cell].append(name)
            else:
                if name not in dict.get(cell):
                    dict[cell].append(name)
                    #dict stores the cell number and then all the reads that had gene assignments
    dictlength=int(cell) #dictionary has as many entries as the last cell counted that had gene assignments (some cells don't have gene assignments)
    for i in range(dictlength+1):
        i=str(i)
        if i not in dict:
            dict[i]=[] #for cells without gene assignments make those cells empty lists in the dict
    for i in range(len(dict)):
        readlist=[]
        i=str(i)
        barcode=bcdict[i] #attach barcode info to cell
        final = open(path + '/assignedreads' + '/'+'cell_'+ i + '_'+barcode+'.new.fasta', 'w') #open file with cell number and barcode
        file=input_path+'/'+'cell_'+i+'_'+barcode+'.fasta'
        seqdict={} 
        for line in open(file): #takes the reads from each cell fasta
            if line.startswith('>'):
                #line=line.rstrip()
                rootname=line.split('_')[0][1:]
                fullname=line.rstrip()[1:] #fetches the full name for each read in these files (so we can determine best consensus read for polishing)
                seqdict[rootname]=[]
                fullnamedict[rootname]=fullname #fullnamedict is cumulative across all single cell fastas 
            else:
                sequence=line.rstrip()
                seqdict[rootname]=sequence #fetches sequence for each read in this file
        for x in dict[i]:
            x=x.split('_')[0]
            name=fullnamedict[x]
            seq=seqdict[x]
            ass=ass_dict[x]
            final.write('%s\t%s\t%s\n' % (name,ass,seq))
        final.close()
    for file in os.listdir(path+'/assignedreads'):
        fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])): #creates cell_#_BC.fasta.new.sorted file with reads sorted by gene_id
        with open(path+'/assignedreads/'+file)as open_file:
            sorted_file=open(path + '/assignedreads/'+ file +'.sorted', 'w')
            rows=open_file.readlines()
            sorted_rows = sorted(rows, key=lambda x: int(x.split()[1][7:]), reverse=False)
            for row in sorted_rows:
                sorted_file.write(row)

def main(input_path, path, ref_genome, gtf_file):
    catfile=concat_fasta(input_path, path)
    #sys.stderr.write('Trimming and concatenating fastas')
    samfile=map_reads(catfile, ref_genome, path)
    sam_con(samfile)
    bamfile=path+'/'+'allreads.bam' 
    fcount(bamfile, gtf)
    countfile=os.getcwd()+'/'+'allreads.bam.featureCounts'
    group_reads(countfile,input_path, path)


main(input_path, path, ref_genome, gtf)
