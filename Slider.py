# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:34:56 2016

@author: daniellefever
"""

#!/usr/bin/python

from __future__ import division



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
