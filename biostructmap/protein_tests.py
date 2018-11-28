"""A collection of tools to handle tests over protein multiple sequence alignments.

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
from .gentests import _calculate_shannon_entropy


def shannon_entropy_protein(alignment, protein_letters=IUPACData.protein_letters,
                            normalized=False, gap='-'):
    '''
    Calculate mean Shannon entropy for all residues in a protein alignment.

    Args:
        alignment (str/Bio.Align.MultipleSequenceAlignment): A multiple sequence
            alignment string in FASTA format or a multiple sequence alignment
            object as a Bio.Align.MultipleSequenceAlignment.
        protein_letters (str, optional): String of all protein letters being
            used to define the amino acid alphabet. Defaults to standard 20
            amino acids. If another alphabet is used (if your sequence contains
            non-standard amino acid), then the maximum Shannon entropy values
            will change accordingly.
        normalized (bool): Normalize such the entropy is in the range [0, 1]

    Returns:
        float: Shannon entropy value.
    '''
    if isinstance(alignment, str):
        alignment = AlignIO.read(StringIO(alignment), format='fasta')
    if not alignment or len(alignment[0]) == 0:
        return None
    translated_positions = list(zip(*[str(x.seq) for
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
