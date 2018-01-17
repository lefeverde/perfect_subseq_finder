# perfect_subseq_finder
Code used to enumerate repetitive sequences from the human genome.

# Identifying Perfect Repetitive Subsequences Within the Human Genome

All scripting was done using python. First code chunk is the imports. The code uses biopython, so that will need to be installed.

Also, the word "perfect" here refers subsequences that are perfect matches.  



```python
from __future__ import division

import string
import numpy 
import os

from Bio import SeqIO, Seq, SeqRecord
from collections import OrderedDict
from collections import Counter


```

The class below is basically just breaks a sequence into subsequences using an adjustable size subsequence size. I don't think it is important to understand the code here. I just used this to break the genome into 20 base pair chunks, for each chromosome. I have included a demonstration of how it works below.


```python
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

```

The code chunk below is just a demonstration of how the sliding window works.


```python
slider_test_str = string.ascii_uppercase
test_slide = SlidingWindow(slider_test_str, 5)
print('Example Sequence:\n')
print(slider_test_str)
print('\nSliding window of subseqs:\n')
for cur_ind,i in enumerate(test_slide):
    pos = cur_ind #-1
    front_astr = '*'*pos
    end_astr = '*'*(len(slider_test_str) - pos - 5)
    out_str = front_astr + i + end_astr
    print out_str
```

    Example Sequence:
    
    ABCDEFGHIJKLMNOPQRSTUVWXYZ
    
    Sliding window of subseqs:
    
    ABCDE*********************
    *BCDEF********************
    **CDEFG*******************
    ***DEFGH******************
    ****EFGHI*****************
    *****FGHIJ****************
    ******GHIJK***************
    *******HIJKL**************
    ********IJKLM*************
    *********JKLMN************
    **********KLMNO***********
    ***********LMNOP**********
    ************MNOPQ*********
    *************NOPQR********
    **************OPQRS*******
    ***************PQRST******
    ****************QRSTU*****
    *****************RSTUV****
    ******************STUVW***
    *******************TUVWX**
    ********************UVWXY*
    *********************VWXYZ


The function below uses the ``SlidingWindow`` to break the genome into 20 bp chunks and then count the number of times each chunk is present. The ambiguous bases are ignored and soft-masked bases are converted to upper case. The subseq argument can be adjusted to an arbitrary length. The ``counter_to_tuplist`` function returns a list with the repetitive sequences sorted by count.



```python
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
```

Here is a code chunk with an example sequence. The following sequence contains the the subsequence 'AATCA' 3 times. The rest of the sequence are bases chosen at random. The 


```python
example_seq = 'AATCAACTAGTTATTACATAATCACAATCA'
counter_to_tuplist(perfect_matcher(example_seq, 5))
```




    [('AATCA', 3),
     ('TACAT', 1),
     ('CAACT', 1),
     ('AGTTA', 1),
     ('CAATC', 1),
     ('TAGTT', 1),
     ('CTAGT', 1),
     ('TAATC', 1),
     ('CATAA', 1),
     ('TTATT', 1),
     ('CACAA', 1),
     ('TTACA', 1),
     ('ACATA', 1),
     ('ATTAC', 1),
     ('AACTA', 1),
     ('ACAAT', 1),
     ('TCAAC', 1),
     ('GTTAT', 1),
     ('TCACA', 1),
     ('TATTA', 1),
     ('ACTAG', 1),
     ('ATAAT', 1),
     ('ATCAA', 1),
     ('ATCAC', 1)]



That is basically the entire script. The remaining function `main` in the `perfect_subseq_finder.py` script takes care of executing the script itself. The script was written to be executed from the command line. The input is a directory containing fasta files of the human genome. I had initially used hg38.p7 but any version would work. Once the files are downloaded and un-compressed, the script can be ran like so: 
```bash
$ python perfect_subseq_finder.py -i fasta_dir
```
