#!/usr/bin/env python3
'''Little script that calculates avg insert length per cell'''

import statistics
import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm
from tqdm import tqdm

parser = argparse.ArgumentParser()
#parser.add_argument('-p', '--input_path', type=str) #path to individual cell .final.fasta UMI-merged/polished files
parser.add_argument('-i', '--input_file', type=str) #where you want 'cellsams' directory to go
#parser.add_argument('-c', '--config_file', type=str) #config file containing paths to executables if not in $PATH
#parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for minimap2 alignment (same used for Map10xUMIs.py)


args = parser.parse_args()
input_file=args.input_file
#input_path=args.input_path
#path=args.output_path
#config_file=args.config_file
#ref_genome=args.ref_genome

'''def avg_length(path):
    fileList=[]
    avg_list=[]
    for file in os.listdir(path):
        fileList.append(file)
    for item in sorted(fileList, key=lambda x: int(x.split('_')[0])):
        value=calc_length(path+'/'+item)
        avg_list.append(value)
    #list_sorted=sorted(avg_list, key=lambda x: int(-x))
    #final_longest=list_sorted[0]
    total_average=statistics.mean(avg_list)
    print(total_average)'''

def calc_length(file):
    len_list=[]
    for line in open(file):
        if line.startswith('@'):
            continue
        else:
            seq=line.rstrip().split('\t')[9]
            print(len(seq))
            #length=len(line)-28
            #len_list.append(length)
    #avg=statistics.mean(len_list)
    #avg_sorted=sorted(len_list, key=lambda x: int(-x))
    #longest_cell=avg_sorted[0]
    #print(avg_sorted)
    #return len_list


#avg_length(input_path)
calc_length(input_file)
#for item in len_list:
    #print(item)
