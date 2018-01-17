from __future__ import division

import argparse
import sys
import os
import subprocess

from Slider import SlidingWindow
from Bio import SeqIO, Seq, SeqRecord
from collections import OrderedDict
from collections import Counter


# Script was taken from:
# https://scipher.wordpress.com/2010/12/02/simple-sliding-window-iterator-in-python/
# This is my original solution.  It is a true iterator class.  It is included here only for completeness. I have replaced it in the repo with the much clearer code above.
class SlidingWindow(object):
    """Returns iterator that will emit chunks of size 'winSize' each time self.next()
    is called."""
    def __init__(self,sequence,winSize,step=1):
        """Returns iterator that will emit chunks of size 'winSize' and 'step' forward in
        the seq each time self.next() is called."""

        # verification code
        if not type(sequence) == type(''):
            raise Exception("**ERROR** type(sequence) must be str.")
        if not ((type(winSize) == type(0)) and (type(step) == type(0))):
            raise Exception("**ERROR** type(winSize) and type(step) must be int.")
        if step > winSize:
            raise Exception("**ERROR** step must not be larger than winSize.")
        if winSize > len(sequence):
            raise Exception("**ERROR** winSize must not be larger than sequence length.")
        self._seq = sequence
        self._step = step
        self._start = 0
        self._stop = winSize

    def __iter__(self):
        return self

    def next(self):
        """Returns next window chunk or ends iteration if the sequence has ended."""
        try:
            assert self._stop <= len(self._seq), "Not True!"
            chunk = self._seq[self._start:self._stop]
            self._start += self._step
            self._stop  += self._step
            return chunk
        except AssertionError:
            raise StopIteration


def perfect_matcher(cur_fasta, subseq_len):
    '''
    This function takes either a single chromosome fasta file
    or a sequence as a string as input. The latter is done for
    testing purposes. If the input (cur_fasta) is not a file, it is assumed
    a string.
    TODO fix error with not existant file as input
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
    opt.add_argument('--no_sorting',
        action='store_true',
        help='turns off sorting of subseqs by count')
    args = parser.parse_args()

    # loop to get all the fastas
    file_list = []
    for cur_file in os.listdir(args.i):
        file_ending = cur_file.split('.')[-1]
        if file_ending == 'fa' or file_ending == 'fasta':
            file_fp = os.path.join(args.i, cur_file)
            file_list.append(file_fp)
    
     # List of input fasta files

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
            if args.no_sorting:
                subseq_cnts = perfect_matcher(cur_fp, args.n)
                # Turns counter to tup list so out objects are the same
                subcounts = [(i[0],i[1]) for i in subseq_cnts.iteritems()]

            else:
                subcounts = counter_to_tuplist(perfect_matcher(cur_fp, args.n))

            #for subseq_str, num in subseq_cnts.iteritems():
            for subseq_str, num in subcounts:
                if num > 1: # dont write singleton subseqs
                    f.write(subseq_str + ',' + str(num) + '\n')

if __name__ == '__main__':
    main()
