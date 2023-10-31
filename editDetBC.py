#!/usr/bin/env python3
'''Trim down detBarcodes.py out file to just BCs that have one or more reads attributed to them'''

import sys
import argparse
#import editdistance as ld #this is for L distance
from tqdm import tqdm
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str)
#parser.add_argument('-n', '--tenx_index_file', type=str)
parser.add_argument('-o', '--output_path', type=str)

args = parser.parse_args()
input_file = args.input_file
#bc_file = args.tenx_index_file
path = args.output_path

def trim_file(input_file, path):
    name='2211_native_detBarcodes_out_edit.fasta'
    out=open(path + name,'w+')
    dict={}
    linenum=1
    for line in open(input_file):
        if line.startswith('>'):
            number=int(line.rstrip().split('_')[2])
            if number > 0:
                barcode=line.rstrip()
                dict[barcode]=[]
        else:
            dict[barcode]=line.rstrip()
        linenum+=1
    for barcode in dict:
        out.write('%s\n%s\n' %(barcode, dict[barcode]))

trim_file(input_file, path)

