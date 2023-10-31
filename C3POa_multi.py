#!/usr/bin/env python3

import argparse
#import editdistance
#import numpy as np
import os
import sys
#import mappy as mm
#import shutil
#from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--splint', type=str) #path to 10x_UMI_splint.fasta
parser.add_argument('-p', '--input_path', type=str) #path to nested fastq directories (i.e /FC1_files)
#parser.add_argument('-o', '--output_path', type=str) #where you want temp output files and assignedreads directory to go
#parser.add_argument('-g', '--gtf_file', type=str) #path to reference genome gtf file for FeatureCount gene id assignments
#parser.add_argument('-c', '--config_file', type=str) #path to config file for minimap2, racon, FeatureCount, and samtools if not in $PATH
#parser.add_argument('-r', '--ref_genome', type=str) #reference GENOME for splice-aware minimap2 alignment




args = parser.parse_args()
splint=args.splint
#input_file=args.input_file
input_path=args.input_path
#path=args.output_path
#config_file=args.config_file
#ref_genome=args.ref_genome
#gtf=args.gtf_file

#C3poa takes the following args ./C3POa.py -r fastq_pass.fastq -o /home/sherbear/ -s 10x_UMI_splint.fasta

list=os.listdir(input_path)
for file in list:
    os.mkdir(input_path+file+'_/')
    reads=input_path+file
    output=input_path+file+'_/'
    os.system('./C3POa-new.py -r %s -o %s -s %s' %(reads, output,splint))

