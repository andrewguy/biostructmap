"""A collection of tools to handle genomic tests.
Part of the structmap package.
"""
import doctest
import dendropy
from structmap.seqtools import (_var_site, _join_alignments,
                                _sliding_window_var_sites)

def tajimas_d(alignment, window=None, step=3):
    '''
    Use DendroPy package to calculate Tajimas D.
    Several optimisations performed to speed up the calculation
    Input needs to be a string representing multiple sequence alignments in
    fasta format.
    Output is Tajima's D value.
    '''
    if window:
        results = []
        prev_win = None
        prev_d = None
        slide = _sliding_window_var_sites(alignment, window, step=step,
                                          isfile=False)
        for i, win in enumerate(slide):
            centre = i*step + 1 + (window-1)/2
            if win == prev_win:
                results.append((centre, prev_d))
            else:
                current_d = _tajimas_d(win)
                results.append((centre, current_d))
                prev_d = current_d
                prev_win = win
        return results
    else:
        #Generate alignment containing only polymorphic sites
        alignment = _join_alignments(_var_site(alignment))
        return _tajimas_d(alignment)

def _tajimas_d(alignment):
    '''
    Use DendroPy to calculate tajimas D.
    '''
    if not alignment or len(alignment[0]) == 0:
        return None
    try:
        seq = dendropy.DnaCharacterMatrix.get(data=alignment.format('fasta'),
                                              schema='fasta')
        taj_d = dendropy.calculate.popgenstat.tajimas_d(seq)
    except TypeError:
        taj_d = None
    return taj_d

if __name__ == '__main__':
    #Test docstrings
    doctest.testmod()
