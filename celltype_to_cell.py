#!/usr/bin/env python3
'''Takes a specified sample ID and WhichCells.txt output from R/Seurat, gets barcodes and fetches matching cell.fasta file from demuxxed dir'''

import argparse
import numpy as np
import os
import sys
import mappy as mm
import shutil
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str)
parser.add_argument('-d', '--demuxed_dir', type=str)
parser.add_argument('-o', '--path', type=str)

args = parser.parse_args()
input_file=args.input_file
demuxed_dir=args.demuxed_dir
path=args.path

'''Takes WhichCells.txt output file from Seurat and fetches barcodes'''
def get_assigned_barcodes(input_file):
    barcode_list=[]
    new_barcode_list=[]
    for line in open(input_file):
        line=line.rstrip().split(' ')
        for i in range(1,len(line)):
            if '-1' in line[i]:
                barcode_list.append(line[i])
    for item in barcode_list:
        identity=item[1:4]
        item=item[5:21]
        if identity == 'NC5':
            new_barcode_list.append(item)
    return new_barcode_list

def get_dir_barcodes(demuxed_dir):
    fileList=[]
    new_fileList=[]
    for file in os.listdir(demuxed_dir):
        if 'final' in file:
            fileList.append(file)
    fileList=sorted(fileList, key=lambda x:(x))
    for item in fileList:
        item=(item.split('.')[0]).split('_')[3]
        new_fileList.append(item)
    #print(new_fileList[8], fileList[8])
    return fileList,new_fileList

def sort_barcodes(assigned_bcs, all_demuxed_bcs, fileList, demuxed_dir,path):
    for i in range (0,len(all_demuxed_bcs)):
        if all_demuxed_bcs[i] in assigned_bcs:
            shutil.copyfile(demuxed_dir+'/'+fileList[i], path+'/'+fileList[i])

assigned_bcs=get_assigned_barcodes(input_file)
fileList, all_demuxed_bcs=get_dir_barcodes(demuxed_dir)
sort_barcodes(assigned_bcs, all_demuxed_bcs, fileList, demuxed_dir,path)

#for item in assigned_bcs:
    #print(item)
