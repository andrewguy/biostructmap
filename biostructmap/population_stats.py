'''A faster implementation than DendroPy.

These methods do not deal with missing or uncertain sequence data.
It is up to the user to filter out uncertain base calls, otherwise methods
will revert to Dendropy functions (which are not guaranteed to be correct for
uncertain base calls).
'''
import numpy as np
import math
from collections import Counter
import operator as op
from functools import reduce

def n_choose_k(n, k):
    "Binomial Coefficient."
    k = min(k, n-k)
    if k < 0:
        return 0
    numer = reduce(op.mul, range(n, n-k, -1), 1)
    denom = reduce(op.mul, range(1, k+1), 1)
    return numer // denom

# Scipy may not be installed. If not, use our own combinatorial function.
# Scipy implementation is slightly faster.
try:
    from scipy.special import comb
except ImportError:
    comb = n_choose_k


def count_differences(seqs):
    '''Calculate the average number of pairwise differences.

    Rather than iterate over all possible pairs of sequences,
    calculates the binomial coefficient for each sequence letter.

    The number of differences is the total number of possible comparisons
    minus the number of comparisons where any given nucleotide letter is
    compared to the same letter.

    Args:
        seqs (list): A list of nucleotide sequences.
    Returns:
        float, int: Average number of pairwise difference, number of
            segregating sites.
    '''
    total_sum = 0
    depth = len(seqs)
    num_seg_sites = 0
    for pos in range(len(seqs[0])):
        base_at_pos = [seq[pos] for seq in seqs]
        count_ = Counter(base_at_pos)
        num_pairwise_differences = comb(depth, 2) - np.sum([comb(count_[base], 2)
                                                            for base in count_])
        total_sum += num_pairwise_differences
        if len(count_) > 1:
            num_seg_sites += 1
    avg_pairwise_diff = total_sum / comb(depth, 2)
    return avg_pairwise_diff, num_seg_sites


def _tajimas_d(num_seqs, avg_pairwise_diffs, num_seg_sites):
    '''Tajima's D formula.

    Args:
        num_seqs (int): The number of sequences.
        avg_pairwise_diffs (float): The average number of pairwise differences.
    Returns:
        float: Tajima's D value.
    '''
    a1 = np.sum(np.reciprocal(np.arange(1, num_seqs, dtype='float64')))
    a2 = np.sum(np.reciprocal(np.square(np.arange(1, num_seqs, dtype='float64'))))
    b1 = float(num_seqs + 1) / (3 * (num_seqs - 1))
    b2 = float(2 * ((num_seqs**2) + num_seqs + 3 )) / (9*num_seqs*(num_seqs-1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - float(num_seqs+2)/(a1 * num_seqs) + float(a2)/(a1 ** 2)
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    d = (
        float(avg_pairwise_diffs - (float(num_seg_sites)/a1))
        / math.sqrt(
            (e1 * num_seg_sites )
          + ((e2 * num_seg_sites) * (num_seg_sites - 1) ))
        )
    return d

def calculate_nucleotide_diversity(seqs):
    avg_pairwise_diff, _num_seg_sites = count_differences(seqs)
    num_sites = len(seqs[0])
    return avg_pairwise_diff / num_sites


def calculate_tajimas_d(seqs):
    num_seqs = len(seqs)
    avg_pairwise_diff, num_seg_sites = count_differences(seqs)
    return _tajimas_d(num_seqs, avg_pairwise_diff, num_seg_sites)

def calculate_wattersons_theta(seqs):
    num_seqs = len(seqs)
    _avg_pairwise_diff, num_seg_sites = count_differences(seqs)
    a1 = np.sum(np.reciprocal(np.arange(1, num_seqs, dtype='float64')))
    return float(num_seg_sites) / a1
