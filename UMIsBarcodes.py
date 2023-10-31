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
parser.add_argument('-t', '--targeted', type=int) #expected cells from 10x experiment
#parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH
#parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for splice-aware minimap2 alignment


args = parser.parse_args()
input_path=args.input_path
targeted=args.targeted
path=args.output_path

'''Makes rank:cell_barcode ordered dictionary by UMI count'''
def UMI_count(input_path):

    fileList=[] 
    for file in os.listdir(input_path):
        #print(file) 
        #file=list(file.split('/'))[len(file)-1]
        if '.final.fasta' in file :
            fileList.append(file)
    #print(fileList)
    #count=0
    cell_list=[]
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):
        count=0
        cell=file.split('_')[1]
        barcode=(file.split('.')[0]).split('_')[2]
        for line in open(input_path+file):
            if line.startswith('>'):
                count +=1
        cellkey=cell+'_'+str(count)+'_'+barcode
        cell_list.append(cellkey) #cell key has cell number_UMIcount_barcode

    rank=0
    sorted_dict={}
    for item in sorted(cell_list, key=lambda x: -int((x.split('_')[1]))): #reorders cellkeys by UMIcount
        number=item.split('_')[1]
        rank+=1
        sorted_dict[rank]=[]
        cell, barcode=item.split('_')[0], item.split('_')[2]
        sorted_dict[rank].append(cell+'_'+barcode)
        sorted_dict[rank].append(number)
        #print('%s\t%s' %(rank, item))
    return sorted_dict

'''Finds UMIcount value of the barcode at 99th percentile and lists all barcodes to be included as cells''' 
def percentile(dict, targeted):
    inclusion_list=[]
    place=targeted-(int(.99*targeted)) #calculates which rank position is at the 99th percentile
    stop=0
    plus=int(targeted*.01) #takes extra %10 of cells
    value=dict[place][1] #fetches the UMI count for the 99th percentile barcode
    inc_value=int(int(value)/10)
    print(value)
    for entry in dict: 
        if int(dict[entry][1]) <= inc_value:
            stop=entry
            break
    fullstop=int(entry+plus) #total number of top-ranking cells to take; any cells exceeding 99th percentile UMI count/10 plus an extra 10% of targeted cells
    for i in range(1,fullstop):
        inclusion_list.append(dict[i][0])

    return inclusion_list


'''Re-sort .final.fasta files by rank for Seurat formatting'''
def resort(list, input_path, path):
    i=0
    zeroes=(len(str(int(len(list)/10))))
    #zeroes=int(len(list)/10)
    
    for item in list:
        file='cell_'+item+'.final.fasta'

        
        if i==len(list):
            break
        elif i < 10:
            shutil.copyfile(input_path+file,path+('0'*zeroes)+str(i)+'_'+file)
        elif i >=10 and i < 100:
            zeroes=(len(str(int(len(list)/10))))-1
            shutil.copyfile(input_path+file,path+('0'*zeroes)+str(i)+'_'+file)
        elif i >= 100 and i < 1000:
            zeroes=(len(str(int(len(list)/10))))-2
            shutil.copyfile(input_path+file,path+('0'*zeroes)+str(i)+'_'+file)
        elif i >= 1000 and i < 10000:
            zeroes=(len(str(int(len(list)/10))))-3
            shutil.copyfile(input_path+file,path+('0'*zeroes)+str(i)+'_'+file)
        elif i >= 10000 and i < 100000:
            zeroes=(len(str(int(len(list)/10))))-4
            shutil.copyfile(input_path+file,path+('0'*zeroes)+str(i)+'_'+file)
        i += 1
    print(len(list))




sorted_dict=UMI_count(input_path)
inclusion_list=percentile(sorted_dict, targeted)
resort(inclusion_list, input_path, path)
