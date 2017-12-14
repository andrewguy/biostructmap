"""A collection of tools to handle genomic tests, in particular Tajima's D.

Part of the biostructmap package.
"""
from __future__ import absolute_import, division, print_function

from io import StringIO
import warnings
from math import log
from Bio import AlignIO
from Bio.Data import IUPACData
from numpy import mean
import dendropy

from .seqtools import _sliding_window_var_sites

def shannon_entropy(alignment, table='Standard',
                    protein_letters=IUPACData.protein_letters,
                    normalized=False, gap='-'):
    '''
    Calculate mean Shannon entropy for all residues in a genomic alignment.

    Args:
        alignment (str/Bio.Align.MultipleSequenceAlignment): A multiple sequence
            alignment string in FASTA format or a multiple sequence alignment
            object, either as a Bio.Align.MultipleSequenceAlignment or a
            biostructmap.SequenceAlignment object.
        table: A codon lookup table used by the Bio.Seq.translate() method.
            See BioPython docs for possible options.
        protein_letters (str, optional): String of all protein letters being
            used to define the amino acid alphabet. Defaults to standard 20
            amino acids. If another alphabet is used (if your sequence contains
            non-standard amino acid), then the maximum Shannon entropy values
            will change accordingly.
        normalized (bool): Normalize such the entropy is in the range [0, 1]
    '''
    if isinstance(alignment, str):
        alignment = AlignIO.read(StringIO(alignment), format='fasta')
    if not alignment or len(alignment[0]) == 0:
        return None
    translated_positions = list(zip(*[str(x.seq.translate(table=table, gap=gap)) for
                                      x in alignment]))
    not_in_alphabet = set([res for x in translated_positions
                           for res in x]).difference(protein_letters)
    if not_in_alphabet:
        warnings.warn("Multiple sequence alignment contains residues that aren't "\
                      "in the provided alphabet. Entropy values will not be "\
                      "accurate - consider supplying an extended amino acid "\
                      "alphabet to the `protein_letters` keyword argument. "\
                      "Offending residue(s) are: {res}".format(res=str(not_in_alphabet)))
    entropy_values = [_calculate_shannon_entropy(x, protein_letters, normalized)
                      for x in translated_positions]
    entropy = mean(entropy_values)
    return entropy

def _calculate_shannon_entropy(seq, protein_letters=IUPACData.protein_letters,
                               normalized=False):
    '''
    Calculate Shannon entropy from a set of residues from a single position.

    Args:
        seq (str): Set of amino acid residues, assuming standard extended alphabet.
        protein_letters (str, optional): String of all protein letters being
            used to define the amino acid alphabet. Defaults to standard 20
            amino acids. If another alphabet is used (if your sequence contains
            non-standard amino acid), then the maximum Shannon entropy values
            will change accordingly.
        normalized (bool): Normalize such the entropy is in the range [0, 1]
    Returns:
        entropy (float): Shannon entropy.
    '''
    # Normalization involves dividing entropy values by the maximum entropy,
    # which is mathematically the same as changing the base of the logarithm
    # to be the size of the protein alphabet.
    if normalized:
        base = len(protein_letters)
    else:
        base = 2
    letters = protein_letters
    entropy = -sum((0 if seq.count(letter) == 0 else
                    seq.count(letter) / len(seq) *
                    log(seq.count(letter)/len(seq), base) for letter in letters))
    return entropy

def tajimas_d(alignment, window=None, step=3):
    """
    Uses DendroPy package to calculate Tajimas D.

    Notes:
        Several optimisations are used to speed up the calculation, including
        memoisation of previous window result, which is used if no new
        polymorphic sites are added/subtracted.

    Args:
        alignment (str/Bio.Align.MultipleSequenceAlignment): A multiple sequence
            alignment string in FASTA format or a multiple sequence alignment
            object, either as a Bio.Align.MultipleSequenceAlignment or a
            biostructmap.SequenceAlignment object.
        window (int, optional): The size of the sliding window over which
            Tajima's D is calculated. Default is None, in which case a
            single Tajima's D value is calculated for the multiple sequence
            alignment.
        step (int, optional): Step size for sliding window calculation.
            Default step size of 3 (ie. one codon).
    Returns:
        float/dict: If window parameter is None, returns a single value for
            Tajima's D. Otherwise a dict mapping genome window midpoint to
            calculated Tajima's D values is returned.
    """
    if window:
        if isinstance(alignment, str):
            alignment = AlignIO.read(StringIO(alignment), 'fasta')
        results = {}
        prev_win = None
        prev_d = None
        slide = _sliding_window_var_sites(alignment, window, step=step)
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
    Uses DendroPy to calculate tajimas D.

    If Tajima's D is undefined (ie. Dendropy Tajima's D method raises a
    ZeroDivisionError), then this method returns None.

    Args:
        alignment (str/Bio.Align.MultipleSequenceAlignment): A multiple sequence
            alignment string in FASTA format or a multiple sequence alignment
            object, either as a Bio.Align.MultipleSequenceAlignment or a
            biostructmap.SequenceAlignment object.

    Returns:
        float: Tajima's D value. Returns None if Tajima's D is undefined.
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

def nucleotide_diversity(alignment):
    """
    Use DendroPy to calculate nucleotide diversity.

    If nucleotide diversity is undefined, returns None.

    Args:
        alignment (str/Bio.Align.MultipleSequenceAlignment): A multiple sequence
            alignment string in FASTA format or a multiple sequence alignment
            object, either as a Bio.Align.MultipleSequenceAlignment or a
            biostructmap.SequenceAlignment object.

    Returns:
        float: Nucleotide diversity value. Returns None if nucleotide
            diversity is undefined.
    """
    if not isinstance(alignment, str):
        data = alignment.format('fasta')
    else:
        data = alignment
    if not alignment or len(alignment[0]) == 0:
        return None
    seq = dendropy.DnaCharacterMatrix.get(data=data, schema='fasta')
    diversity = dendropy.calculate.popgenstat.nucleotide_diversity(seq)
    return diversity

def wattersons_theta(alignment):
    """
    Use DendroPy to calculate Watterson's Theta.

    If Watterson's Theta is undefined, returns None.

    Args:
        alignment (str/Bio.Align.MultipleSequenceAlignment): A multiple sequence
            alignment string in FASTA format or a multiple sequence alignment
            object, either as a Bio.Align.MultipleSequenceAlignment or a
            biostructmap.SequenceAlignment object.

    Returns:
        float: Watterson's Theta value. Returns None if Watterson's Theta is
            undefined.
    """
    if not isinstance(alignment, str):
        data = alignment.format('fasta')
    else:
        data = alignment
    if not alignment or len(alignment[0]) == 0:
        return None
    seq = dendropy.DnaCharacterMatrix.get(data=data, schema='fasta')
    theta = dendropy.calculate.popgenstat.wattersons_theta(seq)
    return theta
