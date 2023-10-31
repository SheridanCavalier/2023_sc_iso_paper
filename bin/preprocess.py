#!/usr/bin/env python3
# Roger Volden

'''Basically we get an output of a readname:splint dict (adapter_dict), the list of splints in the lib(adapter_set), and the number of no splint reads'''
import os
import sys
import mappy as mm
from tqdm import tqdm
import multiprocessing as mp
import shutil
from glob import glob

def preprocess(blat, args, tmp_dir, tmp_adapter_dict, num_reads): #reads that meet the raw read length cutoff in the initial step of C3POa.py are sent here for preprocessing
    tmp_fasta = tmp_dir + 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'
    cat_fasta = tmp_dir + 'concat_blat.fasta'

    # skip the alignment if the psl file already exists
    if not os.path.exists(align_psl) or os.stat(align_psl).st_size == 0: #if the splint blat alignment file does not yet exist, that needs to happen first
        print('Aligning splints to reads with blat', file=sys.stderr)
        chunk_process(num_reads, args, blat) #reads get split into chunks and processed; generating a psl file of splint_to_read_alignments 
    else:
        print('Reading existing psl file', file=sys.stderr)

    #process and chunk_process have been performed on reads, resulting in one large concatenated psl file (align_psl) and fasta file
    
    '''for each read in fasta file, list out all of the alignments for that read from align_psl'''
    blat_dict = {} #dictionary stores all forward and reverse alignments for each read from psl file
    adapter_set = set()
    length = 0
    for line in open(cat_fasta):
        length += 1
    iterator = 0
    infile = open(cat_fasta, 'r')
    while iterator < length:
        line = infile.readline()
        sequence = infile.readline()
        name = line[1:].strip()
        blat_dict[name] = {}
        blat_dict[name]['+'] = []
        blat_dict[name]['-'] = []
        blat_dict[name]['+'].append(('-', 1, 0, 'o'))
        blat_dict[name]['-'].append(('-', 1, len(sequence), 'o'))
        iterator += 2
    #print(blat_dict)
    
    #burn_dict = {}
    for line in open(align_psl):
        a = line.strip().split('\t')
        read_name, adapter, strand = a[9], a[13], a[8]
        gaps, score = float(a[5]), float(a[0])
        sequence_length = int(a[10])
        if gaps < 50 and score > 50: # Looks for unspliced quality alignment
            if strand == '+':
                start = int(a[11]) - int(a[15])
                end = int(a[12]) + int(a[14]) - int(a[16])
            if strand == '-':
                start = int(a[10]) - int(a[12]) #qSize-qEnd
                end = int(a[10]) - int(a[11]) #qSize-qStart
            position = min(max(0, int(start+((end-start)/2))), sequence_length-1)
            blat_dict[read_name][strand].append((adapter, float(a[0]), position, strand)) #blat dict stores all forward and reverse alignment scores and positions for each read
            adapter_set.add(adapter)
    read_dict={}
    for read in mm.fastx_read(args.reads, read_comment=False): #read in original fastq
        readname, seq, qual = read[0], read[1], read[2]
        if readname in blat_dict.keys(): #if the name of the read is in blat_dict
            read_dict[readname]=[]
            read_dict[readname].append(seq)
            read_dict[readname].append(qual)
     
    read_dict, no_splint_reads = create_SW_dict(blat_dict, read_dict, adapter_set, args)
    return read_dict, adapter_set, no_splint_reads
    
def create_SW_dict(blat_dict, read_dict, adapter_set, args):
    SW_dict={}
    no_splint_reads=0
    success_adapter = {}
    for adapter in adapter_set:
        if not os.path.exists(args.out_path + '/' + adapter):
            os.system('mkdir ' + args.out_path + '/' + adapter)
        success_adapter[adapter] = 0
    for read in read_dict:
        name, sequence, quality  = read, read_dict[read][0], read_dict[read][1]
        #print(name, sequence)
        adapter_plus = sorted(blat_dict[name]['+'],
                              key=lambda x: x[1], reverse=True)
        adapter_minus = sorted(blat_dict[name]['-'],        #for each strand, sort alignment entries by score
                             key=lambda x: x[1], reverse=True)

        adapter_all = sorted(adapter_plus + adapter_minus,
                             key=lambda x: x[1], reverse=True) #for each read, sorts alignment entries by score
        #print(name,adapter_plus,adapter_minus)
        best_adapter_direction = adapter_all[0][3] #records the strand for the best scoring alignment for each read
        if best_adapter_direction == 'o':
            no_splint_reads +=1
        plus_list_name, plus_list_position = [], []
        minus_list_name, minus_list_position = [], []
#if there are no alignments for any particular strand, the dict is populated with ('-',1,len,'o')

        for adapter in adapter_plus:
            if adapter[0] != '-': #if plus strand alignment exists...
                plus_list_name.append(adapter[0]) #converts adapter_plus dict into a list with names and positions
                plus_list_position.append(adapter[2])
        for adapter in adapter_minus:
            if adapter[0] != '-': #if minus strand alignment exists...
                minus_list_name.append(adapter[0]) #converts adapter_minus dict into a list with names and positions
                minus_list_position.append(adapter[2])

        plus = False
        if len(plus_list_name) > 0 or len(minus_list_name) > 0: #if there is any splint alignment
            if best_adapter_direction == '+': #and the plus strand is the highest scored alignment...
                plus = True
            if plus:
                adapter = plus_list_name[0]
                SW_dict[name]=[]
                SW_dict[name].append(sequence)
                SW_dict[name].append(quality)
                SW_dict[name].append('+')
                SW_dict[name].append(str(plus_list_position[0]))
            else:
                adapter = minus_list_name[0] #the adapter is going to be '10x_UMI_Splint' for any of our reads with alignments
                SW_dict[name]=[]
                SW_dict[name].append(sequence)
                SW_dict[name].append(quality)
                SW_dict[name].append('-')
                SW_dict[name].append(str(minus_list_position[0]))
            splint_reads_folder = adapter
            success_adapter[adapter] += 1 #counts how many reads align to each splint (for us again it's just 10x_UMI_Splint)
            success = success_adapter[adapter] #if the best splint alignment was to the reverse strand, make note of that for revcomp SW fxn
    return SW_dict, no_splint_reads
    
    
def cat_files(path, pattern, output): #we may be able to use this to also concat the tmp_fastas 
    '''Use glob to get around bash argument list limitations'''
    final_psl = open(output, 'w+')
    for f in tqdm(glob(path + pattern), desc='Catting files'):
        with open(f) as psl:
            for line in psl:
                final_psl.write(line)
    final_psl.close()

    
def remove_files(path, pattern):
    '''Use glob to get around bash argument list limitations'''
    for d in tqdm(glob(path + pattern), desc='Removing preprocessing files'):
        shutil.rmtree(d)

        
def process(args, reads, blat, iteration): #called in the chunk_process fxn with args=(args, tmp_reads, blat, iteration)
    tmp_dir = args.out_path + 'pre_tmp_' + str(iteration) + '/' #makes tmp dir for this particular chunk of reads
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    tmp_fa = tmp_dir + 'tmp_for_blat.fasta'
    tmp_fa_fh = open(tmp_fa, 'w+') #write temp fasta file
    for header, seq in reads.items(): #for readname:seq key:entry pair in the tmp_reads dictionary...
        print('>' + header, file=tmp_fa_fh) #generate temp fasta file from this tmp_reads name:seq dictionary for blat alignment
        print(seq, file=tmp_fa_fh)
    tmp_fa_fh.close()
    align_psl = tmp_dir + 'tmp_splint_aln.psl' #create psl from this chunk
    b_msgs = tmp_dir + 'blat_messages.log'

    os.system('{blat} -noHead -stepSize=1 -t=DNA -q=DNA -minScore=15 -minIdentity=10 {splint} {reads} {psl} >{blat_msgs}'
              .format(blat=blat, splint=args.splint_file, reads=tmp_fa, psl=align_psl, blat_msgs=b_msgs)) #perform blat alignment using splint argument -s against reads in this chunk
    #os.remove(tmp_fa)
    #the psl file for this particular chunk is now stored in 'tmp_splint_aln.psl'
    
    
def chunk_process(num_reads, args, blat): #chunks input and processes
    '''Split the input fasta into chunks and process'''
    if args.blatThreads: #takes argument -b put into C3POa which determines number of blat threads to use for multiprocessing
        chunk_size = (num_reads // args.numThreads) + 1 #divide reads into number of threads
    else:
        chunk_size = args.groupSize #otherwise the default group size is used for each chunk
    if chunk_size > num_reads:
        chunk_size = num_reads

    pool = mp.Pool(args.numThreads)
    pbar = tqdm(total=num_reads // chunk_size + 1, desc='Preprocessing')
    iteration, current_num, tmp_reads, target = 1, 0, {}, chunk_size #set tmp_reads as empty dict
    for read in mm.fastx_read(args.reads, read_comment=False): #args.reads is our raw fastq_pass infile that we put into C3POa
        if len(read[1]) < args.lencutoff: #filter out reads that don't make length cutoff -l 
            continue
        tmp_reads[read[0]] = read[1] #populate tmp_read dictionary at key=read name with entry=read sequence
        current_num += 1
        if current_num == target:
            pool.apply_async(process, args=(args, tmp_reads, blat, iteration), callback=lambda _: pbar.update(1)) #perform process fxn on this chunk
            iteration += 1
            target = chunk_size * iteration
            if target >= num_reads:
                target = num_reads
            tmp_reads = {}
    pool.close()
    pool.join()
    pbar.close()

    cat_files(args.out_path,'pre_tmp_*/tmp_splint_aln.psl',args.out_path + 'tmp/splint_to_read_alignments.psl') #all of the tmp_splint_aln.psl files from the read chunks are concatted together into the final splint_to_read_alignments.psl file
    cat_files(args.out_path, 'pre_tmp_*/tmp_for_blat.fasta', args.out_path + 'tmp/concat_blat.fasta') #all passing reads catted into fasta file
    remove_files(args.out_path, 'pre_tmp*')
