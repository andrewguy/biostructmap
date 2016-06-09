"""A collection of tools to handle genomic tests.
Part of the structmap package.
"""
from __future__ import absolute_import, division, print_function

import dendropy
from io import StringIO
from Bio import AlignIO
from .seqtools import (_var_site, _join_alignments,
                       _sliding_window_var_sites)

def tajimas_d(alignment, window=None, step=3):
    """
    Use DendroPy package to calculate Tajimas D.
    Several optimisations performed to speed up the calculation
    Input needs to be a string representing multiple sequence alignments in
    fasta format or a SeqIO alignment object.
    Output is Tajima's D value.
    """
    if window:
        if isinstance(alignment, str):
            alignment = AlignIO.read(StringIO(alignment), 'fasta')
        results = {}
        prev_win = None
        prev_d = None
        slide = _sliding_window_var_sites(alignment, window, step=step,
                                          isfile=False)
        for i, win in enumerate(slide):
            centre = i*step + 1 + (window-1)/2
            if win == prev_win:
                results[centre] = prev_d
            else:
                current_d = _tajimas_d(win)
                results[centre] = current_d
                prev_d = current_d
                prev_win = win
        return results
    else:
        return _tajimas_d(alignment)

def _tajimas_d(alignment):
    """
    Use DendroPy to calculate tajimas D.
    """
    if not isinstance(alignment, str):
        data = alignment.format('fasta')
    else:
        data = alignment
    if not alignment or len(alignment[0]) == 0:
        return None
    try:
        seq = dendropy.DnaCharacterMatrix.get(data=data,
                                              schema='fasta')
        taj_d = dendropy.calculate.popgenstat.tajimas_d(seq)
    except ZeroDivisionError:
        taj_d = None
    return taj_d
