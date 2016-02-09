"""Helper module that contains various utility functions for structmap
package.
"""
from __future__ import absolute_import, division, print_function

import re

def to_string(sequence):
    '''
    Takes a Bio.Seq object and converts to a sequence string.
    '''
    return ''.join([x for x in sequence])


def is_unambig_dna(text, search=re.compile(r'[^ATGCatgc]').search):
    '''
    Checks if sequence string is unambiguous DNA.
    Returns false if the string contains any character other than ATGC.
    Usage:
    >>> is_unambig_dna("ACTGTCTGTCAT")
    True
    >>> is_unambig_dna("ACTGTCGCMCTCG")
    False
    '''
    return not bool(search(text))
