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
parser.add_argument('-p', '--input_path', type=str) #path to directory with re-ordered/ranked final.fasta files
parser.add_argument('-b', '--barcode_guide', type=str) #path to bcGuide.txt created by demux_nano.py script
parser.add_argument('-o', '--output_path', type=str) #expected cells from 10x experiment
#parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH
#parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for splice-aware minimap2 alignment

args = parser.parse_args()
input_path=args.input_path
barcode_guide=args.barcode_guide
path=args.output_path

filelist=[]
for file in os.listdir(input_path):
    filelist.append(file)

guide=path+'/BcGuide_reordered.txt'
write_guide=open(guide, 'w')
for file in sorted(filelist, key=lambda x: int(x.split('_')[0])):
    barcode=(file.split('_')[3]).split('.')[0]
    for line in open(barcode_guide):
        if barcode in line:
            cell,bar=line.split('\t')[0], line.split('\t')[1]
            write_guide.write(cell+'\t'+bar)

