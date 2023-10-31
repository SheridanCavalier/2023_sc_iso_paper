#!/usr/bin/env python3
'''Wrapper for multiple'''
'''Sorting fastq reads into demuxed cell files to create cell_01.fastq etc...'''


import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-q', '--subreads', type=str) #path to dir with subread fastqs
parser.add_argument('-p', '--input_path', type=str) #path to demuxed dir

args = parser.parse_args()
subreads=args.subreads
input_path=args.input_path

list=os.listdir(subreads)
for file in list:
    if 'fastq' in file:
        subs=subreads+file
        output=input_path
        os.system('./make_cell_subreads.py %s %s %s' %(input_path,subs,output))

