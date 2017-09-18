'''A collection of tools to handle sequence manipulation.
Part of the biostructmap package.
'''
from __future__ import absolute_import, division, print_function

from io import StringIO
import operator
import re
import subprocess
import tempfile
from Bio import AlignIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.pairwise2 import align
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#Use local BLAST+ installation. Falls back to pairwise2 if False.
LOCAL_BLAST = True
#Use local exonerate installation to align dna to protein sequences.
#Falls back to a basic method using either BLAST+ or pairwise2 if False,
#but won't take into consideration introns or frameshift mutations.
LOCAL_EXONERATE = True

def _sliding_window(seq_align, window, step=3, fasta_out=False):
    '''
    Generate a Multiple Sequence Alignment over a sliding window.

    Input is either filehandle, or Bio.AlignIO multiple sequence alignment
    object.
    Args:
        seq_align: A multiple sequence alignment. Either a filehandle, or
            Bio.AlignIO multiple sequence alignment object.
        window (int): Sliding window width
        step (int, optional): Step size to increment each window. Default of 3.
        fasta_out (bool): If True, output will be a fasta formatted string. If
            False, then output will be an AlignIO object.
    Yields:
        str/MultipleSequenceAlignment: The next window in the sliding window
            series for the original multiple sequence alignment.
    '''
    try:
        alignments = AlignIO.read(seq_align, 'fasta')
    except AttributeError:
        alignments = seq_align
    #Length of alignments
    length = len(alignments[0])
    for i in range(0, length-window, step):
        alignment = alignments[:, i:i+window]
        if fasta_out:
            alignment = alignment.format('fasta')
        yield alignment

def _sliding_window_var_sites(seq_align, window, step=3):
    '''
    Generate a Multiple Sequence Alignment over a sliding window, only
    including polymorphic sites in the alignment.

    Notes:
        Returns an empty MultipleSequenceAlignment object if no polymorphic
        sites are found within the window.

    Args:
        seq_align: A multiple sequence alignment. Either a filehandle, or
            Bio.AlignIO multiple sequence alignment object.
        window (int): Sliding window width
        step (int, optional): Step size to increment each window. Default of 3.
        fasta_out (bool): If True, output will be a fasta formatted string. If
            False, then output will be an AlignIO object.
    Yields:
        MultipleSequenceAlignment: The next window in the sliding window
            series for the original multiple sequence alignment, with
            only polymorphic sites displayed.
    '''
    try:
        alignments = AlignIO.read(seq_align, 'fasta')
    except AttributeError:
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
    # Add/remove sites from the end/start of window as appropriate.
    for i in range(0, (length-window), step):
        for j in range(step):
            if i + j in align_dict:
                alignment = alignment[:, 1:]
            if i + j + window in align_dict:
                alignment = alignment + align_dict[i+j+window]
        yield alignment


def _var_site(alignment):
    '''
    Take a multiple sequence alignment object and return polymorphic sites in a
    dictionary object.

    This function is used to simplify the input to a tajima's D calculation.

    Args:
        alignment: A multiple sequence alignment object.

    Returns:
        dict: A dictionary containing polymorphic sites (value) accessed by
            position in the alignment (key).
    '''
    result = {}
    for i in range(len(alignment[0])):
        site = alignment[:, i]
        #Check if string contains a single character. Most efficient method
        #found so far.
        if site != len(site) * site[0]:
            result[i] = alignment[:, i:i+1]
    return result

def _join_alignments(align_dict):
    '''
    Take a dictionary of multiple sequence alignments, and join according to
    dictionary key order (generally position in sequence).

    Args:
        align_dict (dict): A dictionary containing single-site multiple sequence
            alignment objects accessed by position in original alignment.

    Returns:
        MultipleSequenceAlignment: A multiple sequence alignment object
            containing all polymorphic sites.
    '''
    output = None
    for key in sorted(align_dict):
        if not output:
            output = align_dict[key]
        else:
            output = output + align_dict[key]
    return output

def align_protein_sequences(comp_seq, ref_seq):
    '''
    Perform a pairwise alignment of two sequences.

    Uses BLAST+ if LOCAL_BLAST is set to True, otherwise uses Bio.pairwise2.

    Args:
        comp_seq (str): A comparison protein sequence.
        ref_seq (str): A reference protein sequence.

    Returns:
        dict: A dictionary mapping comparison sequence numbering (key) to
            reference sequence numbering (value)
        dict: A dictionary mapping reference sequence numbering (key) to
            comparison sequence numbering (value)
    '''
    if LOCAL_BLAST:
        return blast_sequences(comp_seq, ref_seq)
    else:
        return pairwise_align(comp_seq, ref_seq)


def align_protein_to_dna(prot_seq, dna_seq):
    '''
    Aligns a protein sequence to a genomic sequence.

    If LOCAL_EXONERATE flag is set to True, takes into consideration
    introns, frameshifts and reverse-sense translation if using Exonerate.

    If LOCAL_EXONERATE flag is set to False, then a simple translation and
    pairwise alignment is performed, and does not consider introns, frameshifts
    or reverse-sense translations.

    Note:
        This method uses the external program Exonerate:
        http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
        This needs to be installed in the users PATH.

    Args:
        prot_seq (str): A protein sequence.
        dna_seq (str): A genomic or coding DNA sequence

    Returns:
        dict: A dictionary mapping protein residue numbers to codon positions:
            e.g. {3:(6,7,8), 4:(9,10,11), ...}
    '''
    if LOCAL_EXONERATE:
        return _align_prot_to_dna_exonerate(prot_seq, dna_seq)
    else:
        return _align_prot_to_dna_no_exonerate(prot_seq, dna_seq)



def _align_prot_to_dna_no_exonerate(prot_seq, dna_seq):
    '''
    Aligns a protein sequence to a genomic sequence. Does not take consider
    introns, frameshifts or reverse-sense translation.
    If these are required, should use Exonerate method instead.

    Args:
        prot_seq (str): A protein sequence.
        dna_seq (str): A genomic or coding DNA sequence

    Returns:
        dict: A dictionary mapping protein residue numbers to codon positions:
            e.g. {3:(6,7,8), 4:(9,10,11), ...}
    '''
    #Translate DNA sequence to protein sequence
    dna_prot_seq = str(Seq(dna_seq, generic_dna).translate())
    #Use existing methods to align protein-protein
    prot_dna_dict, _ = align_protein_sequences(prot_seq, dna_prot_seq)
    #Convert output to protein: codon dict
    protein_to_codons = {key: (value*3-2, value*3-1, value*3) for
                         key, value in prot_dna_dict.items()}
    return protein_to_codons

def pairwise_align(comp_seq, ref_seq):
    '''
    Perform a pairwise alignment of two sequences.

    Uses the BioPython pairwise2 module with the BLOSUM62 matrix for scoring
    similarity. Gap opening penalty is -11 and gap extend penalty is -1,
    which is the same as the default blastp parameters.

    Output is two dictionaries: residue numbering in PDB chain (key) mapped to
    the residue position in the reference sequence (value), and vice versa.

    Args:
        comp_seq (str): A comparison protein sequence.
        ref_seq (str): A reference protein sequence.

    Returns:
        dict: A dictionary mapping comparison sequence numbering (key) to
            reference sequence numbering (value)
        dict: A dictionary mapping reference sequence numbering (key) to
            comparison sequence numbering (value)
    '''
    alignment = align.globalds(comp_seq, ref_seq, matlist.blosum62, -11, -1,
                               penalize_end_gaps=False, one_alignment_only=True)[0]
    query_string = alignment[0]
    sbjct_string = alignment[1]
    #Create dictionary mapping position in PDB chain to position in ref sequence
    pdb_to_ref = {}
    ref_to_pdb = {}
    key = 1
    ref = 1
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

def blast_sequences(comp_seq, ref_seq):
    '''
    Perform BLAST of two protein sequences using NCBI BLAST+ package.

    Output is two dictionaries: residue numbering in PDB chain (key) mapped to
    the residue position in the reference sequence (value), and vice versa.

    Notes:
        User must have NCBI BLAST+ package installed in user's PATH.

    Args:
        comp_seq (str): A comparison protein sequence.
        ref_seq (str): A reference protein sequence.

    Returns:
        dict: A dictionary mapping comparison sequence numbering (key) to
            reference sequence numbering (value)
        dict: A dictionary mapping reference sequence numbering (key) to
            comparison sequence numbering (value)
    '''
    with tempfile.NamedTemporaryFile(mode='w') as comp_seq_file, \
         tempfile.NamedTemporaryFile(mode='w') as ref_seq_file:
        comp_seq_file.write(">\n" + str(comp_seq) + "\n")
        ref_seq_file.write(">\n" + str(ref_seq) + "\n")
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


def _construct_sub_align_from_chains(alignments, codons, fasta=False):
    '''
    Take a list of biostructmap multiple sequence alignment objects, and
    return a subset of codons based on an input list in the form
    [('A',(1,2,3)),('B',(4,5,6)),...].

    Notes:
        Codons should be 1-indexed, not 0-indexed.

    Args:
        alignment (dict): A dictionary of multiple sequence alignment objects
            accessed by a tuple of chain ids for each alignment.
        codons (list): a subset of codons in a list of the form [(1,2,3),...]
        fasta (bool, optional): If True, will return multiple sequence
            alignment as a string in FASTA format.

    Returns:
        MulitpleSequenceAlignment: A subset of the initial alignment as a
            multiple sequence alignment object. If the fasta kwarg is set to
            True, returns a string instead.
    '''
    chain_alignments = {}
    chain_strains = {}
    for key, alignment in alignments.items():
        chain_alignments[key] = alignment.get_alignment_position_dict()
        chain_strains[key] = alignment.get_isolate_ids()
    codons = [(chain_id, x) for chain_id, sublist in codons for x in sublist]
    sub_align = []
    for codon in codons:
        #List is zero indexed, hence the need to call codon-1
        sub_align.append(list(chain_alignments[codon[0]][codon[1]-1]))
    _sub_align_transpose = zip(*sub_align)
    sub_align_transpose = [''.join(x) for x in _sub_align_transpose]
    if fasta:
        strains = list(chain_strains.values())[0]
        fasta_out = ''.join('>{}\n{}\n'.format(*t) for t in
                            zip(strains, sub_align_transpose))
        return fasta_out
    return sub_align_transpose

def _construct_sub_align(alignment, codons, fasta=False):
    '''
    Take a biostructmap multiple sequence alignment object, and return a
    subset of codons based on an input list in the form [(1,2,3),(4,5,6),...].

    Notes:
        Codons should be 1-indexed, not 0-indexed.

    Args:
        alignment: A multiple sequence alignment object.
        codons (list): a subset of codons in a list of the form [(1,2,3),...]
        fasta (bool, optional): If True, will return multiple sequence
            alignment as a string in FASTA format.

    Returns:
        MulitpleSequenceAlignment: A subset of the initial alignment as a
            multiple sequence alignment object. If the fasta kwarg is set to
            True, returns a string instead.
    '''
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
        fasta_out = ''.join('>{}\n{}\n'.format(*t) for t in
                            zip(strains, sub_align_transpose))
        return fasta_out
    return sub_align_transpose


def _align_prot_to_dna_exonerate(prot_seq, dna_seq):
    '''
    Aligns a protein sequence to a genomic sequence. Takes into consideration
    introns, frameshifts and reverse-sense translation.

    Note:
        This method uses the external program Exonerate:
        http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
        This needs to be installed in the users PATH.

    Args:
        prot_seq (str): A protein sequence.
        dna_seq (str): A genomic or coding DNA sequence

    Returns:
        dict: A dictionary mapping protein residue numbers to codon positions:
            e.g. {3:(6,7,8), 4:(9,10,11), ...}
    '''
    #TODO Use Biopython exonerate parser. Didn't realise that existed when I wrote this parser.
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
        raise UserWarning("Exonerate direction doesn't match either '+' or '-'")

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
            for _ in range(count2):
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
            for _ in range(count2):
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
