"""A collection of tools to handle sequence manipulation.
Part of the structmap package.
"""
from __future__ import absolute_import, division, print_function

from Bio import AlignIO, pairwise2 as pw2
from Bio.SubsMat import MatrixInfo

def _sliding_window(seq_align, window, step=3, isfile=True, fasta_out=False):
    """
    Generator function that takes a multiple sequence alignment, and generates a
    Multiple Sequence Alignment over a sliding window.
    Input is either filehandle, or Bio.AlignIO multiple sequence alignment
    object.
    Set isfile=False if dealing with Bio.AlignIO MSA object.
    Output is either AlignIO object or fasta string.
    If fasta_out==False, then output will be AlignIO object.
    """
    if isfile:
        alignments = AlignIO.read(seq_align, 'fasta')
    else:
        alignments = seq_align
    #Length of alignments
    length = len(alignments[0])
    for i in range(0, length-window, step):
        alignment = alignments[:, i:i+window]
        if fasta_out:
            alignment = alignment.format('fasta')
        yield alignment

def _sliding_window_var_sites(seq_align, window, step=3, isfile=True):
    """
    Generator function that takes a multiple sequence alignment,
    and generates a Multiple Sequence Alignment over a sliding window, only
    including polymorphic sites in the alignment.
    Output is Bio.AlignIO alignment object.
    """
    if isfile:
        alignments = AlignIO.read(seq_align, 'fasta')
    else:
        alignments = seq_align
    #Length of alignments
    length = len(alignments[0])

    align_dict = _var_site(alignments)

    #Create first window
    initial_sites = {key:value for (key, value) in align_dict.items()
                     if key < window}
    #Small hack to set type of 'initial_sites' variable if no alignments fall
    #within initial window
    initial_sites[-1] = alignments[:, 0:0]

    alignment = _join_alignments(initial_sites)
    yield alignment
    #Add/remove sites from the end/start of window as appropriate.
    for i in range(0, (length-window), step):
        for j in range(step):
            if i + j in align_dict:
                alignment = alignment[:, 1:]
            if i + j + window in align_dict:
                alignment = alignment + align_dict[i+j+window]
        yield alignment


def _var_site(alignment):
    """
    Take a multiple sequence alignment object and return polymorphic sites in a
    dictionary object.
    Use this function to simplify the input to a tajima's D calculation.
    """
    result = {}
    for i in range(len(alignment[0])):
        site = alignment[:, i]
        #Check if string contains a single character. Most efficient method
        #found so far.
        if site != len(site) * site[0]:
            result[i] = alignment[:, i:i+1]
    return result

def _join_alignments(align_dict):
    """
    Take a dictionary of multiple sequence alignments, and join according to
    dictionary key order (generally position in sequence).
    """
    output = None
    for key in sorted(align_dict):
        if not output:
            output = align_dict[key]
        else:
            output = output + align_dict[key]
    return output

def map_to_sequence(compseq, refseq):
    """
    Takes a protein refence sequence and corresponding PDB file sequence as
    input, and maps the residue numbering of the PDB chain to the reference
    sequence, using the Bio.pairwise2 package for alignment.
    Output is two dictionaries: residue numbering in PDB chain (key) mapped to
    the residue position in the reference sequence (value), and vice versa.
    """
    matrix = MatrixInfo.blosum62
    gap_open = -10
    gap_extend = -0.5

    if not isinstance(compseq, str):
        compseq = ''.join([x for x in compseq])
    if not isinstance(refseq, str):
        refseq = ''.join([x for x in refseq])

    alns = pw2.align.globalds(compseq, refseq, matrix, gap_open, gap_extend)
    aln_key, aln_ref, _, _, _ = alns[0]
    #Create dictionary mapping position in PDB chain to position in ref sequence
    pdb_to_ref = {}
    ref_to_pdb = {}
    key = 0
    ref = 0
    for i, res in enumerate(aln_key):
        if res.isalpha() and aln_ref[i].isalpha():
            key += 1
            ref += 1
            pdb_to_ref[key] = ref
            ref_to_pdb[ref] = key
        elif res.isalpha():
            key += 1
        elif aln_ref[i].isalpha():
            ref += 1

    return pdb_to_ref, ref_to_pdb

def _construct_sub_align(alignments, codons):
    """
    Take a multiple sequence alignment object, and return a subset of codons
    based on an input list in the form [(1,2,3),(4,5,6),...].
    Return subset of the initial alignment as a multiple sequence alignment
    object.
    """
    codons = [x for sublist in codons for x in sublist]
    sub_align = {}
    #Small hack to set type of 'sub_align' variable if no alignments fall
    #within initial window
    sub_align[-1] = alignments[:, 0:0]
    for codon in codons:
        sub_align[codon] = alignments[:, codon:codon+1]
    output = _join_alignments(sub_align)
    return output

def prot_to_dna_position(dna_indices, prot_indices):
    """
    Take a list of nucleotide positions in a reference sequence, and a map
    of protein residue numbers matching the nucleotide positions,
    and return a dictionary with key: residue number (int), value: nucleotide
    positions (list of 3 values). Eg. {4:[7,8,9]}
    >>> prot_indices = [5,6,8]
    >>> dna_indices = range(1,10)
    >>> x = prot_to_dna_position(dna_indices,prot_indices)
    >>> print([(y,x[y]) for y in sorted(x)])
    [(5, (1, 2, 3)), (6, (4, 5, 6)), (8, (7, 8, 9))]
    """
    lookup_dict = {x:tuple(dna_indices[i*3:(i+1)*3]) for i, x in
                   enumerate(prot_indices)}
    return lookup_dict
