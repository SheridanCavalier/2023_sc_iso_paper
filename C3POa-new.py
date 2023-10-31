#!/usr/bin/env python3
# Roger Volden
# Sheridan Cavalier edits 2021

import os
import sys
import numpy as np
import argparse
import multiprocessing as mp
import mappy as mm
from conk import conk
from tqdm import tqdm
import gc
import gzip
import shutil
from glob import glob

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

from preprocess import preprocess
from call_peaks import call_peaks
from determine_consensus import determine_consensus

VERSION = 'v2.2.3--self to self'

def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Makes consensus sequences from R2C2 reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--reads', '-r', type=str, action='store',
                          help='FASTQ file that contains the long R2C2 reads.')
    parser.add_argument('--splint_file', '-s', type=str, action='store',
                          help='Path to the splint FASTA file.')
    parser.add_argument('--out_path', '-o', type=str, action='store', default=os.getcwd(),
                        help='''Directory where all the files will end up.
                                Defaults to your current directory.''')
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='''If you want to use a config file to specify paths to
                                programs, specify them here. Use for racon and blat
                                if they are not in your path.''')
    parser.add_argument('--lencutoff', '-l', type=int, action='store', default=1000,
                        help='''Sets the length cutoff for your raw sequences. Anything
                                shorter than the cutoff will be excluded. Defaults to 1000.''')
    parser.add_argument('--mdistcutoff', '-d', type=int, action='store', default=500,
                        help='''Sets the median distance cutoff for consensus sequences.
                                Anything shorter will be excluded. Defaults to 500.''')
    parser.add_argument('--zero', '-z', action='store_false', default=True,
                        help='Use to exclude zero repeat reads. Defaults to True (includes zero repeats).')
    parser.add_argument('--numThreads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--groupSize', '-g', type=int, default=1000,
                        help='Number of reads processed by each thread in each iteration. Defaults to 1000.')
    parser.add_argument('--blatThreads', '-b', action='store_true', default=False,
                        help='''Use to chunk blat across the number of threads instead of by groupSize (faster).''')
    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')
    parser.add_argument('--version', '-v', action='version', version=VERSION, help='Prints the C3POa version.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

def configReader(path, configIn):
    progs = {}
    with open(configIn) as f:
        for line in f:
            if line.startswith('#') or not line.rstrip().split():
                continue
            line = line.rstrip().split('\t')
            progs[line[0]] = line[1]
    possible = set(['racon', 'blat'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
    # check for missing programs
    # if missing, default to path
    for missing in possible - inConfig:
        path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

def cat_files(path, pattern, output, description, compress):
    '''Use glob to get around bash argument list limitations'''
    if compress:
        output += '.gz'
        final_fh = gzip.open(output, 'wb+')
    else:
        final_fh = open(output, 'w+')
    for f in tqdm(glob(path + pattern), desc=description):
        with open(f) as fh:
            for line in fh:
                if compress:
                    line = line.encode()
                final_fh.write(line)
def cat_files(path, pattern, output, description):
    '''Use glob to get around bash argument list limitations'''
    final_fh = open(output, 'w+')
    for f in tqdm(glob(path + pattern), desc=description):
        with open(f) as fh:
            for line in fh:
                final_fh.write(line)
    final_fh.close()

def remove_files(path, pattern):
    '''Use glob to get around bash argument list limitations'''
    for d in tqdm(glob(path + pattern), desc='Removing files'):
        shutil.rmtree(d)

def rounding(x, base):
    '''Rounds to the nearest base, we use 50'''
    return int(base * round(float(x) / base))

'''Reverse comp fxn'''
def reverse_complement(sequence): #revcomp for minus strand alignments
  Seq = ''
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
  for item in sequence[::-1]:
    Seq = Seq + complement[item]
  return Seq


def analyze_reads(args, tmp_reads, adapter_set, iteration, racon): #analyze_reads fxn called on reads in chunks (tmp_reads that have met raw length criteria)
    penalty, iters, window, order = 20, 3, 41, 2
    for read in tmp_reads:
        name=read
        seed, seq, strand, qual = int(tmp_reads[read][3]), tmp_reads[read][0], tmp_reads[read][2], tmp_reads[read][1]
        seq_len=len(seq)
        #print (name) #, seed, seq, strand)
        if strand == '-':
            seq=reverse_complement(seq)
            seed=seq[seed:seed+700]
            #splint=reverse_complement(splint)
        else:
            seed=seq[seed:seed+700]
            
        scores = conk.conk(seed, seq, penalty) #feeding conk the seed, the sequence and the penalty to generate scores
        #scores2 = conk.conk(splint, seq, penalty)
        peaks = call_peaks(scores, 50 , iters, window, order) #take 50bases as min subread dist
        #peaks2 = call_peaks(scores2, 500 , iters, window, order)
        if not list(peaks):
            continue

        '''Convert peak positions to subread start positions''' 
        #adj = -(len(seed)//2)+100 #formula is unique to this version of C3POa
        peaks = list(peaks + 100) #convert peaks to approx subread starts

        for i in range(len(peaks) - 1, -1, -1): #count backwards from list of peaks and delete peaks before start of read or after end of read
            if peaks[i] >= seq_len:
                del peaks[i]
        if not peaks:
            continue

        subreads, qual_subreads = [], []
        if len(peaks) > 1:
            subread_lens = np.diff(peaks)
            subread_lens = [rounding(x, 50) for x in subread_lens]
            median_subread_len = np.median(subread_lens)
            for i in range(len(subread_lens)):
                bounds = [peaks[i], peaks[i+1]]
                if median_subread_len*0.8 <= subread_lens[i] <= median_subread_len*1.2:
                    subreads.append(seq[bounds[0]:bounds[1]])
                    qual_subreads.append(qual[bounds[0]:bounds[1]])
            if peaks[0] > 100:
                subreads.append(seq[:peaks[0]])
                qual_subreads.append(qual[:peaks[0]])
            if seq_len - peaks[-1] > 100:
                subreads.append(seq[peaks[-1]:])
                qual_subreads.append(qual[peaks[-1]:])
        else:
            subreads.append(seq[:peaks[0]])
            qual_subreads.append(qual[:peaks[0]])
            subreads.append(seq[peaks[0]:])
            qual_subreads.append(qual[peaks[0]:])

        tmp_dir = args.out_path + '10x_UMI_Splint_tmp' + str(iteration) + '/'
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        subread_file = tmp_dir + 'subreads.fastq'

        consensus, repeats = determine_consensus(
            args, read, subreads, qual_subreads,
            racon, tmp_dir, subread_file
        )


        if consensus:
            avg_qual = round(sum([ord(x)-33 for x in qual])/seq_len, 2)
            cons_len = len(consensus)
            final_out = open(tmp_dir + '/R2C2_Consensus.fasta', 'a+')
            print('>' + name + '_' + '_'.join([str(x) for x in [avg_qual, seq_len, repeats, cons_len]]), file=final_out)
            print(consensus, file=final_out)
            final_out.close()


 
        '''subreads, qual_subreads, dangling_subreads, qual_dangling_subreads = [], [], [], []
        if len(peaks) > 1:
            subread_lens = np.diff(peaks) #returns list of diffs between peaks
            subread_lens = [rounding(x, 50) for x in subread_lens] #round off subreads lens to the nearest multiple of 50
            median_subread_len = np.median(subread_lens) #takes median bc non-normal dist
            #median_subread_len1 = np.median(subread_lens1)
            for i in range(len(subread_lens)):
                bounds = [peaks[i], peaks[i+1]]
                if median_subread_len*0.8 <= subread_lens[i] <= median_subread_len*1.2: #if subread len is within margin of error around median
                    subreads.append(seq[bounds[0]:bounds[1]]) #take sequence of subread
                    qual_subreads.append(qual[bounds[0]:bounds[1]]) #take per-base qual scores of subread
            if peaks[0] > 100: #identifies partial subreads at beginning of raw read
                dangling_subreads.append(seq[:peaks[0]])
                qual_dangling_subreads.append(qual[:peaks[0]])
            if seq_len - peaks[-1] > 100: #identifies partial subreads at end of raw read
                dangling_subreads.append(seq[peaks[-1]:])
                qual_dangling_subreads.append(qual[peaks[-1]:])
        else: #if only one peak is found, identify the partial subreads before and after it
            dangling_subreads.append(seq[:peaks[0]])
            qual_dangling_subreads.append(qual[:peaks[0]])
            dangling_subreads.append(seq[peaks[0]:])
            qual_dangling_subreads.append(qual[peaks[0]:])

        tmp_dir=args.out_path + '10x_UMI_Splint_tmp' + str(iteration) + '/'
        if not os.path.isdir(tmp_dir):
           os.mkdir(tmp_dir)
        subread_file = tmp_dir + 'subreads.fastq'
        
        consensus, repeats = determine_consensus(args,read, subreads, qual_subreads, dangling_subreads, qual_dangling_subreads, racon, tmp_dir, subread_file)


        if consensus:
            avg_qual = round(sum([ord(x)-33 for x in qual])/seq_len, 2)
            cons_len = len(consensus)
            final_out = open(tmp_dir + '/R2C2_Consensus.fasta', 'a+')
            print('>' + name + '_' + '_'.join([str(x) for x in [avg_qual, seq_len, repeats, cons_len]]), file=final_out)
            print(consensus, file=final_out)
            final_out.close()'''

def main(args):
    if not args.out_path.endswith('/'):
        args.out_path += '/'
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)
    log_file = open(args.out_path + 'c3poa.log', 'w+')

    if args.config: #if config file given
        progs = configReader(args.out_path, args.config)
        racon = progs['racon']
        blat = progs['blat']
    else:
        racon = 'racon' #otherwise racon and blat are the respective executables already in $PATH
        blat = 'blat'

    tmp_dir = args.out_path + 'tmp/' #make 'tmp/' directory in outdir
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    # read in the file and preprocess
    read_list, total_reads = [], 0
    short_reads = 0
    tmp_fasta = tmp_dir + 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'

    tmp_adapter_dict = {}
    #reads in raw fastq_pass seq file
    for read in mm.fastx_read(args.reads, read_comment=False): #mappy.fastx_read returns tuple for each read that contains {name,seq,qual}
        if len(read[1]) < args.lencutoff: #separates raw reads that are below raw read cutoff specified by -l (default=1000)
            short_reads += 1 #count reads that don't satisfy length requirement
            continue
        tmp_adapter_dict[read[0]] = [[None, 1, None]] # [adapter, matches, strand]
        total_reads += 1
    read_dict, adapter_set, no_splint_reads = preprocess(blat, args, tmp_dir, tmp_adapter_dict, total_reads) #calls preprocessing fxn from preprocess.py
    #basically we get an output of a readname:splint dict, the list of splints in the lib, and the number of no splint reads
    #preprocess separates reads by the which splint aligns best, this is only useful for splint multiplexing--for us we want preprocessing to yield the best alignment position and strand of our 10x_UMI_splint in each read for seeding in SW

    
    for adapter in adapter_set:
        if not os.path.exists(args.out_path + adapter):
            os.mkdir(args.out_path + adapter) #make separate dirs for demultiplexed splints (we will only have 10x_UMI_Splint)
    all_reads = total_reads + short_reads #caclutes percentage of reads that don't meet raw length cutoff and don't have a splint alignment (via blat)
    print('C3POa version:', VERSION, file=log_file)
    print('Total reads:', all_reads, file=log_file)
    print('No splint reads:',
           no_splint_reads,
           '({:.2f}%)'.format((no_splint_reads/all_reads)*100),
           file=log_file)
    print('Under len cutoff:',
           short_reads,
           '({:.2f}%)'.format((short_reads/all_reads)*100),
           file=log_file)
    print('Total thrown away reads:',
           short_reads + no_splint_reads,
           '({:.2f}%)'.format(((short_reads + no_splint_reads)/all_reads)*100),
           file=log_file)
    print('Reads after preprocessing:', all_reads - (short_reads + no_splint_reads), file=log_file)
    log_file.close()

    

    pool = mp.Pool(args.numThreads, maxtasksperchild=1)
    pbar = tqdm(total=total_reads // args.groupSize + 1, desc='Calling consensi')
    iteration, current_num, tmp_reads, target = 1, 0, {}, args.groupSize
    #print(len(read_dict))
    for read in read_dict: #mapp.fastx_read (fastq, False) returns a tuple for each read that stores {read name, seq, qual}
        tmp_reads[read]=[]
        tmp_reads[read].append(read_dict[read][0])
        tmp_reads[read].append(read_dict[read][1])
        tmp_reads[read].append(read_dict[read][2])
        tmp_reads[read].append(read_dict[read][3])
    #print(tmp_reads)
    #print(len(tmp_reads))
    #analyze_reads(args, tmp_reads, adapter_set, iteration, racon)
        current_num +=1 #count each read in the dictionary  
        if current_num == target:
            pool.apply_async(analyze_reads,args=(args, tmp_reads, adapter_set, iteration, racon), callback=lambda _: pbar.update(1)) 
            iteration += 1
            target = args.groupSize * iteration
            if target >= total_reads:
                target = total_reads
            tmp_reads = {}
            gc.collect()
    pool.close()
    pool.join()
    pbar.close()

    
    cat_files(args.out_path, '/10x_UMI_Splint_tmp*/R2C2_Consensus.fasta',args.out_path + '/R2C2_Consensus.fasta','Catting consensus reads')
    cat_files(args.out_path, '/10x_UMI_Splint_tmp*/subreads.fastq',args.out_path + '/R2C2_Subreads.fastq','Catting subreads')
    remove_files(args.out_path, '/10x_UMI_Splint_tmp*')

if __name__ == '__main__':
    args = parse_args()
    if not args.reads or not args.splint_file:
        print('Reads (--reads/-r) and splint (--splint_file/-s) are required', file=sys.stderr)
        sys.exit(1)
    mp.set_start_method("spawn")
    main(args)

