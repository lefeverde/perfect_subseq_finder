from __future__ import division

import argparse
import sys
import os
import subprocess

from Slider import SlidingWindow
from Bio import SeqIO, Seq, SeqRecord
from collections import OrderedDict
from collections import Counter


def perfect_matcher(cur_fasta, subseq_len):
    '''
    Counts the number of times a given subsequence
    occurs in a sequence
    '''

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
    '''
    Function which takes care of executing the script. I wrote
    this script to be executed from the command line. Hence, all
    the parser stuff below.

    The script is meant to be run from the command line like so:
    $ python perfect_subseq_finder.py -i fasta_dir

    The -n argument specificies the subsequence length. The default is
    set to 20.

    '''
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


    file_list = os.listdir(args.i) # List of input fasta files

    for cur_file in file_list:
        # Loop which iterates over files
        cur_fp = os.path.join(args.i, cur_file)
        # automagically generates output file name
        out_handle = cur_file.split('.')[0] + '_perfect_subseq_cnts.csv'
        # prints status since this script can take a while
        sys.stdout.write('Working on file: ' + cur_file + '\n')
        sys.stdout.flush()
        # opens output file for writing
        with open(out_handle, 'w+') as f:

            subseq_cnts = perfect_matcher(cur_fp, args.n)
            for subseq_str, num in subseq_cnts.iteritems():
                f.write(subseq_str + ',' + str(num) + '\n')


if __name__ == '__main__':
    main()
