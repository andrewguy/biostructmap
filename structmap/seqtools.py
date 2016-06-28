"""A collection of tools to handle sequence manipulation.
Part of the structmap package.
"""
from __future__ import absolute_import, division, print_function

import subprocess
import tempfile
import re
import operator
from io import StringIO
from Bio import AlignIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

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


def blast_sequences(comp_seq, ref_seq):
    '''
    Perform BLAST of two protein sequences using NCBI BLAST+ package.
    Output is two dictionaries: residue numbering in PDB chain (key) mapped to
    the residue position in the reference sequence (value), and vice versa.
    '''
    if not isinstance(comp_seq, str):
        comp_seq = ''.join([x for x in comp_seq])
    if not isinstance(ref_seq, str):
        ref_seq = ''.join([x for x in ref_seq])

    with tempfile.NamedTemporaryFile(mode='w') as comp_seq_file, \
         tempfile.NamedTemporaryFile(mode='w') as ref_seq_file:
        comp_seq_file.write(">\n" + comp_seq + "\n")
        ref_seq_file.write(">\n" + ref_seq + "\n")
        ref_seq_file.flush()
        comp_seq_file.flush()
        blastp_cline = NcbiblastpCommandline(query=comp_seq_file.name,
                                             subject=ref_seq_file.name,
                                             evalue=0.001, outfmt=5)
        alignment, _stderror = blastp_cline()
    blast_xml = StringIO(alignment)
    blast_record = NCBIXML.read(blast_xml)
    temp_score = 0
    high_scoring_hsp = None
    #Retrieve highest scoring HSP
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.score > temp_score:
                temp_score = hsp.score
                high_scoring_hsp = hsp
    query_string = high_scoring_hsp.query
    sbjct_string = high_scoring_hsp.sbjct
    #Create dictionary mapping position in PDB chain to position in ref sequence
    pdb_to_ref = {}
    ref_to_pdb = {}
    key = high_scoring_hsp.query_start
    ref = high_scoring_hsp.sbjct_start
    for i, res in enumerate(query_string):
        if res.isalpha() and sbjct_string[i].isalpha():
            pdb_to_ref[key] = ref
            ref_to_pdb[ref] = key
            key += 1
            ref += 1
        elif res.isalpha():
            key += 1
        elif sbjct_string[i].isalpha():
            ref += 1
    return pdb_to_ref, ref_to_pdb


def _construct_sub_align(alignment, codons, fasta=False):
    """
    Take a structmap multiple sequence alignment object, and return a
    subset of codons based on an input list in the form [(1,2,3),(4,5,6),...].
    Return subset of the initial alignment as a multiple sequence alignment
    object.
    Note that codons should be 1-indexed.
    """
    alignments = alignment.get_alignment_position_dict()
    strains = alignment.get_isolate_ids()
    codons = [x for sublist in codons for x in sublist]
    sub_align = []
    for codon in codons:
        #List is zero indexed, hence the need to call codon-1
        sub_align.append(list(alignments[codon-1]))
    _sub_align_transpose = zip(*sub_align)
    sub_align_transpose = [''.join(x) for x in _sub_align_transpose]
    if fasta:
        fasta_out = ''.join('>{}\n{}\n'.format(*t) for t in zip(strains, sub_align_transpose))
        return fasta_out
    return sub_align_transpose

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

def align_protein_to_dna(prot_seq, dna_seq):
    """
    Aligns a protein sequence to a genomic sequence. Takes into consideration
    introns frameshifts and reverse-sense translation.

    Note: This method uses the external program Exonerate:
          http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
          This needs to be installed in the users Path.

    Output from align_protein_to_dna is a dictionary mapping protein residue
    numbers to codon positions: {3:(6,7,8), 4:(9,10,11), ...}
    """
    with tempfile.NamedTemporaryFile(mode='w') as protein_seq_file, \
         tempfile.NamedTemporaryFile(mode='w') as dna_seq_file:
        protein_seq_file.write(">\n" + prot_seq + "\n")
        dna_seq_file.write(">\n" + dna_seq + "\n")
        dna_seq_file.flush()
        protein_seq_file.flush()
        #If protein sequence length is small, then exonerate score needs
        #to be adjusted in order to return alignment.
        #With a length n, a perfect match would score 5n.
        #Hence we make a threshold of 3n (60%).
        exonerate_call = ["exonerate",
                          "--model", "protein2genome",
                          "--showalignment", "False",
                          "--showvulgar", "True",
                          protein_seq_file.name,
                          dna_seq_file.name]
        if len(prot_seq) < 25:
            threshold = str(len(prot_seq) * 3)
            exonerate_call.append("--score")
            exonerate_call.append(threshold)
        alignment = subprocess.check_output(exonerate_call)
    vulgar_re = re.search(r"(?<=vulgar:).*(?=\n)",
                              alignment.decode("utf-8"))
    if not vulgar_re:
        raise UserWarning("Did not find exonerate alignment.")
    vulgar_format = vulgar_re.group(0)
    protein_start = vulgar_format.split()[0]
    dna_start = vulgar_format.split()[3]
    matches = vulgar_format.split()[7:]
    direction = vulgar_format.split()[5]
    protein_count = int(protein_start)
    dna_count = int(dna_start)

    if direction == "+":
        step = operator.add
    elif direction == "-":
        step = operator.sub
        dna_count += 1
    else:
        raise UserWarning("Exonerate direction does not match either '+' or '-'")


    if len(matches) % 3:
        raise UserWarning("The vulgar output from exonerate has failed \
                           to parse correctly")
    #Split output into [modifier, query_count, ref_count] triples
    matches = [matches[i*3:i*3+3] for i in range(len(matches)//3)]
    matched_bases = {}

    codon = []

    #Convert vulgar format to dictionary with residue: codon pairs
    for region in matches:
        modifier = region[0]
        count1 = int(region[1])
        count2 = int(region[2])
        if modifier == 'M':
            if count1 != count2 / 3:
                raise UserWarning("Match in vulgar output is possibly " +
                                  "incorrect - number of protein residues " +
                                  "should be the number of bases divided by 3")
            for i in range(count2):
                dna_count = step(dna_count, 1)
                codon.append(dna_count)
                if len(codon) == 3:
                    protein_count += 1
                    matched_bases[protein_count] = tuple(codon)
                    codon = []
        if modifier == 'C':
            if count1 != count2 / 3:
                raise UserWarning("Codon in vulgar output is possibly " +
                                  "incorrect - number of protein residues " +
                                  "should be the number of bases divided by 3")
            raise UserWarning("Unexpected output in vulgar format - not " +
                              "expected to need functionality for 'codon' " +
                              "modifier")
        if modifier == 'G' or modifier == 'N':
            if codon:
                raise UserWarning("Warning - split codon over gap in " +
                                  "exonerate output!")
            protein_count = protein_count + count1
            dna_count = step(dna_count, count2)
        if modifier == '5' or modifier == '3':
            if count1 != 0:
                raise UserWarning("Warning - protein count should be 0 in " +
                                  "exonerate output over intron splice sites.")
            dna_count = step(dna_count, count2)
        if modifier == 'I':
            if count1 != 0:
                raise UserWarning("Warning - protein count should be 0 in " +
                                  "exonerate output over intron.")
            dna_count = step(dna_count, count2)
        if modifier == 'S':
            for i in range(count2):
                dna_count = step(dna_count, 1)
                codon.append(dna_count)
                if len(codon) == 3:
                    protein_count += 1
                    matched_bases[protein_count] = tuple(codon)
                    codon = []
        if modifier == 'F':
            raise UserWarning("Unexpected frameshift in exonerate output - " +
                              "check alignment input.")

    return matched_bases
