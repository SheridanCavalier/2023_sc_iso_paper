#!/usr/bin/env python3

'''creates barcodes.tsv, features.tsv and matrix file for transcript counts--uses salmon quant.sf files for each cell in given directory'''

import os
import sys
import argparse


'''Input args'''
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to directory with individual cell salmon quant.sf files
parser.add_argument('-o', '--output_path', type=str) #where you want the output directory to go
#parser.add_argument('-', '--config_file', type=str) #config file containing paths to executables if not in $PATH
#parser.add_argument('-r', '--ref_tx', type=str) #reference TRANSCRIPTOME for minimap2 alignment


args = parser.parse_args()
input_path=args.input_path
path=args.output_path
#config_file=args.config_file
#ref_tx=args.ref_tx

def parsedir(input_path,outPath):
    '''Makes barcodes.tsv, barcode:position dictionary for all cell files in input_dir and makes transcript list from first salmon.sf file'''
    barcode_dict={}
    fileList=[]
    pos=0
    num_features=117494
    bcOut = open(outPath + 'barcodes.tsv', 'w+')
    mtxOut = open(outPath + 'matrix.mtx', 'w+')
    mtxOut.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    '''After matrix file is made you need to ADD the second line which totals 1) numfeatures 2) numcells 3)num lines in matrix file'''

    for file in os.listdir(input_path):
        fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[0])):
        pos += 1
        bc=file.rstrip().split('_')[3]
        barcode=bc.split('.')[0]
        bcOut.write(barcode + '-1\n')
        barcode_dict[pos]=barcode

        line_num=-1
        feature_dict={}
        if pos==1:
            txDict=maketx(file,input_path,outPath)
        for line in open(input_path+'/'+file):
            line_num += 1
            if line.startswith('Name'):
                continue
            numReads=float(line.rstrip().split('\t')[4])
            if numReads > 0:
                mtxOut.write(str(line_num) + ' ' + str(pos) + ' ' + str(round(numReads)) + '\n')
    mtxOut.close()
    return barcode_dict, txDict



def maketx(inFile,input_path,path): #makes dictionary with transcript:position key:value pairs
    '''
    Makes gene dictionary out of the fist quant.sf file in the input directory and creates 'features.tsv' file with transcriptID. Tx dictionary stores the line number
    in features.tsv where each gene ID is for matrix file formatting.
    '''
    genesOut = open(path + 'features.tsv', 'w+')
    txDict = {}
    pos=1
    for line in open(input_path+'/'+inFile):
        if line.startswith('Name'):
            continue
        tx_id = line.rstrip().split('\t')[0]
        genesOut.write(tx_id + '\t' + tx_id+ '\t'+'Gene Expression'+ '\n')
        txDict[tx_id]=pos
        pos +=1
    genesOut.close()
    return txDict

def makeBC(inFile):
    '''
    Makes the barcodes.tsv file by adding '-1' to the end of each barcode--these numbers change if multiple 10x wells were used!
    '''
    numCells = 0
    bcOut = open(outPath + 'barcodes.tsv', 'w+')
    sys.stderr.write('Writing barcodes to {}barcodes.tsv\n'.format(outPath))
    for line in open(inFile):
        line = line.rstrip().split('\t')
        bcOut.write(line[1] + '-1\n')
        numCells += 1
    bcOut.close()
    return numCells


def modifyFC(fcIn): #this basically gets the counts per gene (we want a way to store the cell number info for each count column too I would think)
    #We can make a second dictionary that holds all the cell names and the index of cell name in this dict vs the gene count id ones can be the same. 
    '''
    Reads the featureCounts output and throws it into a dictionary
    countDict = {geneID: [0, 1, 0, 0, ...], ...}
    The index of the list +1 is the cell # #we can just have the cell number be the column header because that is what it is in mine
    '''
    first = True
    countDict = {}
    sys.stderr.write('Reading in the featureCounts output\n')
    for line in open(fcIn): #for each line in the featurecounts file...
        if line.startswith('#'): # featureCounts command info
            continue
        if first:
            first = False
            cellnames=line.rstrip().split('\t')[6:]
            cellnumindex=[]
            for i in range(0,len(line.rstrip().split('\t')[6:]),1):
                cellnumber=((cellnames[i].split('/')[1]).split('.')[0]).split('_')[0]
                cellnumindex.append(cellnumber) #creates cellnum list that has the cell number in the same index position as the genecount dict entries
            continue
        #print(cellnumindex)
        line = line.rstrip().split('\t')
        gene = line[0]
        counts = line[6:] # only need the count columns
        countDict[gene] = counts #creates dict with gene:count/cell key:value pairs. Entries have same indexing as cellnumindex so counts can be mapped back to the right cell
    return countDict, cellnumindex

def make_mtx(cellnumindex,countDict,geneDict):
    mtxOut = open(outPath + 'matrix.mtx', 'w+')

    genecounts=0 #number of lines in mtx file
    for i in range(0,len(cellnumindex),1):
        for entry in countDict:
            if countDict[entry][i] != '0':
                genecounts += 1
                #print(geneDict[entry],i+1,countDict[entry][i])
    print(len(geneDict)) #number of cells
    mtxOut = open(outPath + 'matrix.mtx', 'w+')
    sys.stderr.write('Writing matrix to {}matrix.mtx\n'.format(outPath))
    mtxOut.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    mtxOut.write("{0} {1} {2}\n".format(str(len(geneDict)), str(len(cellnumindex)), str(genecounts)))

    for i in range(0,len(cellnumindex),1):
        for entry in countDict:
            if countDict[entry][i] != '0':
                mtxOut.write(str(geneDict[entry]) + ' ' + str(i+1) + ' ' + str(countDict[entry][i]) + '\n')
 




barcode_dict, txDict = parsedir(input_path,path)

#for key in barcode_dict:
   
#print(barcode_dict)


#countDict,cellnumindex=modifyFC(exp) #takes counts per gene per cell info from featurecounts file and puts into dictionary
#geneDict = makeGenes(anno) #makes gene dictionary from annotation file and creates features.tsv

#make_mtx(cellnumindex,countDict,geneDict) #uses cell number index, counts per gene, and gene dictionary to make MEX file

