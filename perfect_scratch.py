# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:34:56 2016

@author: daniellefever
"""

#!/usr/bin/python

from __future__ import division
import string
import argparse
import sys
import os
import subprocess

from pathlib2 import Path
from Slider import SlidingWindow
from Bio import SeqIO, Seq, SeqRecord
from collections import OrderedDict
from collections import Counter
from itertools import izip

parser = argparse.ArgumentParser()
req = parser.add_argument_group('required arguments')
opt = parser.add_argument_group('optional arguments')

req.add_argument('-i',
    help = 'dir of fasta files',
    metavar = '',
    required=True)
opt.add_argument('-n',
    type=int,
    default = 20,
    help='sliding window length (Default: 20)',
    metavar='')


args = parser.parse_args()


def perfect_matcher(cur_fasta, subseq_len):
    '''
    This function takes either a single chromosome fasta file
    or a sequence as a string as input. The latter is done for
    testing purposes. If the input (cur_fasta) is not a file, it is assumed
    a string.
    TODO fix error with not existent file as input
    '''

    if os.path.isfile(cur_fasta):
        with open(cur_fasta, 'rU') as open_fasta:
            cur_seqrec = SeqIO.read(open_fasta, 'fasta')
            seq_as_str = cur_seqrec.seq
    else:
        seq_as_str = cur_fasta

    electric_slide = SlidingWindow(str(seq_as_str), subseq_len)
    big_subseq_list = []
    for i in electric_slide:
        if 'N' not in i:
            big_subseq_list.append((i.upper()))

    #big_subseq_list = [i.upper() for i in electric_slide]
    subseq_cnts = Counter(big_subseq_list)
    return subseq_cnts


"""
def perfect_matcher(cur_fasta, subseq_len):


    with open(cur_fasta, 'rU') as open_fasta:
        cur_seqrec = SeqIO.read(open_fasta, 'fasta')
        electric_slide = SlidingWindow(str(cur_seqrec.seq), subseq_len)
        big_subseq_list = []
        for i in electric_slide:
            if 'N' not in i:
                big_subseq_list.append((i.upper()))

        #big_subseq_list = [i.upper() for i in electric_slide]
        subseq_cnts = Counter(big_subseq_list)
    return subseq_cnts
"""


def counter_to_tuplist(counter_object):
    '''
    Takes counter object from perfect_matcher and
    turns it into a sorted list of tuples in the format as
    (subseq, count)
    '''
    subcount_list = []
    for subseq, count in counter_object.iteritems():
        subcount_list.append((subseq, count))
    sorted_subcount = sorted(subcount_list, key=lambda x: -x[1])
    return sorted_subcount



def main():

    file_list = os.listdir(args.i)
    for cur_file in file_list:
        cur_fp = os.path.join(args.i, cur_file)
        out_handle = cur_file.split()[0] + '_perfect_subseq_cnts.csv'
        sys.stdout.write('Working on file: ' + cur_file + '\n')
        sys.stdout.flush()
        with open(out_handle, 'w+') as f:
            subseq_cnts = perfect_matcher(cur_fp, args.n)
            for subseq_str, num in subseq_cnts.iteritems():
                f.write(subseq_str + ',' + str(num) + '\n')


base_chr21 = [i.split(',') for i in open('raw_fastas/chr21.fa_perfect_subseq_cnts.csv')]
del(base_chr21)

base_chr21 = []
for i in  open('raw_fastas/chr21.fa_perfect_subseq_cnts.csv'):
    seq, num = i.split(',')
    num = int(num)
    if num > 10:
        base_chr21.append((seq, str(num)))

base_chr21 = sorted(base_chr21,key=lambda x: -int(x[1]))
base_chr21[0:100]

original_out = '/Users/daniellefever/Desktop/subseq_finder/original_outputs_from_dovarak/chr21.fa_perfect_subseq_cnts.csv'
new_out = '/Users/daniellefever/Desktop/subseq_finder/raw_fastas/chr21.fa_perfect_subseq_cnts.csv'
oo = [i for i in open(original_out)]
no = [i for i in open(new_out)]


fq = '/Users/daniellefever/Downloads/Rotation_Related/Xing_Lab/crisprability/raw_data/fastas/chr21.fa'
cur_seqrec = SeqIO.read(fq, 'fasta')
seq_str = str(cur_seqrec.seq)
len(seq_str)
seq_rm_ambigs = seq_str.replace('N', '')
cs = seq_rm_ambigs[0:500]

electric_slide = SlidingWindow(str(seq_str), 500)


for cur_sub in electric_slide:
    cs = cur_sub #kludge to get around weird str error
    num_ambigs = cs.count('N')
    if num_ambigs < 200 and num_ambigs > 1:
    #if str.islower(cur_sub):

        break


ambigs_subseq="""NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGatcttcctccaaagaaattgtagttttcttctggcttagaggtagatcatcttggtccaatcagactgaaatgccttgaggctagatttcagtctttgtggcagctggtgaatttctagtttgccttttcagctagggattagctttttaggggtcccaatgcctagggagatttctaggtcctctgttccttgctgacctccaattttgtctatccttttgctgagaggtctgcttaacttccttttagtcaggtagctccattttatgctaagcttcttagttgctcaccttctgcag"""

lower_subseq= """tgtggtgttatgatatatattggttttcatccacagttcctggcttataactcccctagcacttgttacagtcttttgttataatattgggtgtattaggcctcaggagcaggcctctcaccttctcatggccttttttcatttttatgttcctgcctttctggttgtgggtcttaagaccatctcaagagagagtcccaccctataccctggagggaggaatgctgatatcatgaaacttccataaaaatccaggaggacagggttcagtgagcttctgggtagttgaacacatggatgttcctgtagggtggcccgcccagggatggcatggaagctctgctcccttcccctatgaattgctctaagtgtccttcatctatatcctttgcaatgtcctttataaaacaccaggaaatgtaagtgtttccctgagttctgtgagccactccaacaaattaatcaaacccaaagagggggtcctgagaagccaactcgaa"""
upper_subseq = """AGGAAAACTTTAAATAGCATTGAGTTATCAGTACTTTCATGTCTTGATACATTTCTTCTTGAAAATGTTCATGCTTGCTGATTTGTCTGTTTGTTGAGAGGAGAATGTTCAGAATTTTATATCTTCAACATCTTTTTCTTCATTAATAAGATACTGAGATTTTATAACTCTTGTCATTTTGGTCACTTATATTTTCATATGGAAATATCGTATAATCCAGGGTTTCCAATATATTTGTGTAAAATTAAGAAAATTATCTTATCTAATAACTTGATCAATATCTGTGATTATATTTTCATTGCCTTCCAATTTTAATATTTGTTCTCTATTCCTTCTTAATCTGGATTGAAGTTCTGATTAATTATTTTAATGTTGCAAATTGTTTTCACTTTTTCCATAAAATGAGTTCTAGAGTTTATTTCTTTACTGCATCATTCTATTTTCAAGTCATGAACTTCTGCTTCAACTAAAAAAAAAAAACTCACCGTTTGTATGAAATT"""





out_handle = os.path.join('test_fastas', 'ambig_subseq.fa')
subseqs = [ambigs_subseq, lower_subseq, upper_subseq]
out_names = ['ambigs_subseq', 'lower_subseq', 'upper_subseq']
for i in izip(out_names, subseqs):
    out_handle = os.path.join('test_fastas', i[0] + '.fa')
    with open(out_handle, 'w+') as f:
        header = '>' + i[0] + '\n'
        f.write(header)
        f.write(i[1])



from_str = counter_to_tuplist(perfect_matcher(ambigs_subseq, 4))
from_file = counter_to_tuplist(perfect_matcher(out_handle, 4))
assert(from_str == from_file)


slider_test_str = string.ascii_uppercase
test_slide = SlidingWindow(slider_test_str, 5)
for cur_ind,i in enumerate(test_slide):
    pos = cur_ind #-1
    front_astr = '*'*pos
    end_astr = '*'*(len(slider_test_str) - pos)
    out_str = front_astr + i + end_astr
    print out_str

if __name__ == '__main__':
    main()
