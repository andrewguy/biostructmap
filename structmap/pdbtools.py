"""A collection of tools for analysing a pdb file.

Helper module for the structmap package.
"""
from __future__ import absolute_import, division, print_function

from Bio.SeqIO import PdbIO
from Bio.SeqUtils import seq1, ProtParamData
from Bio.Data.SCOPData import protein_letters_3to1
from scipy.spatial import distance
import numpy as np
from .seqtools import _construct_sub_align, blast_sequences
from . import gentests

ss_lookup_dict = {
    'H': 0,
    'B': 1,
    'E': 2,
    'G': 3,
    'I': 4,
    'T': 5,
    'S': 6,
    '-': 7,
    0: 'H',
    1: 'B',
    2: 'E',
    3: 'G',
    4: 'I',
    5: 'T',
    6: 'S',
    7: '-'
    }

def _euclidean_distance_matrix(chain, selector='all'):
    """Compute the Euclidean distance matrix for all atoms in a pdb chain.

    Args:
        chain (Chain): Bio.PDB Chain object.
        selector (str): The atom in each residue with which to compute
            distances. The default setting is 'all', which gets all
            non-heterologous atoms. Other potential options include 'CA', 'CB'
            etc. If an atom is not found within a residue object, then method
            reverts to using 'CA'.
    Returns:
        np.array: A euclidean distance matrix.
        np.array: A reference list of all atoms in the model (positionally matched to the euclidean matrix).
    """
    reference = []
    coords = []
    #Filter on non-HET atoms
    for residue in (x for x in chain if x.get_id()[0] == ' '):
        #If selecting based on all atoms within residue
        if selector == 'all':
            for atom in residue:
                coords.append(atom.get_coord())
                reference.append(atom.get_full_id()[3][1])
        #If measuring distance on particular atoms
        else:
            if selector in residue:
                select_atom = selector
            #Revert to carbon alpha if atom is not found
            elif 'CA' in residue:
                select_atom = 'CA'
            #if CA is not found, do not include residue in distance matrix
            else:
                continue
            coords.append(residue[select_atom].get_coord())
            reference.append(residue[select_atom].get_full_id()[3][1])
    #Convert to a np array, and compute Euclidean distance.
    coord_array = np.array(coords)
    euclid_mat = distance.pdist(coord_array, 'euclidean')
    #Convert to squareform matrix
    euclid_mat = distance.squareform(euclid_mat)
    ref_array = np.array(reference)
    return euclid_mat, ref_array

def nearby(chain, radius=15, selector='all'):
    """
    Takes a Bio.PDB chain object, and find all residues within a radius of a
    given residue.

    Args:
        chain (Chain): Bio.PDB Chain object.
        radius (float/int): The radius (Angstrom) over which to select nearby
            residues
        selector (str): The atom in each residue with which to compute
            distances. The default setting is 'all', which gets all
            non-heterologous atoms. Other potential options include 'CA', 'CB'
            etc. If an atom is not found within a residue object, then method
            reverts to using 'CA'.
    Returns:
        dict: A dictionary containing nearby residues for each
            residue in the chain.
    """
    #Setup variables
    ref_dict = {}
    euclidean_distance, reference = _euclidean_distance_matrix(chain, selector)
    within_radius = euclidean_distance <= radius
    near_map = within_radius * reference
    #Iterate over all atoms in Euclidean distance matrix.
    for i, atom in enumerate(near_map):
        if atom[i] not in ref_dict:
            ref_dict[atom[i]] = atom[np.nonzero(atom)]
        else:
            ref_dict[atom[i]] = np.append(ref_dict[atom[i]],
                                          atom[np.nonzero(atom)])
    #Remove non-unique values.
    for key in ref_dict:
        ref_dict[key] = np.unique(ref_dict[key])
    return ref_dict


def get_pdb_seq(filename):
    """
    Get a protein sequence from a PDB file.

    Will return multiple sequences if PDB file contains several chains.

    Args:
        filename (str/filehandle): A PDB filename or file-like object.
    Returns:
        dict: Protein sequences (str) accessed by chain id.
    """
    #Open PDB file and get sequence data
    try:
        with open(filename, 'r') as f:
            seq = [s for s in PdbIO.PdbSeqresIterator(f)]
    except TypeError:
        #If file-like object is passed instead (io.StringIO)
        seq = [s for s in PdbIO.PdbSeqresIterator(filename)]
    #A bit of manipulation to get Seq object into a dictionary
    #Key is chain ID, and value is sequence as a string.
    try:
        sequences = {s.id.split(":")[1]:''.join([x for x in s]) for s in seq}
    except IndexError:
        sequences = {s.id:''.join([x for x in s]) for s in seq}
    return sequences

def get_pdb_seq_from_atom(chain):
    """
    Get a protein sequence from chain atoms in a PDB file.

    This is used as a 'last resort' when sequence is not available in PDB
    headers.

    Args:
        chain: A Bio.PDB chain object.
    Returns:
        str: Protein sequence.
    """
    seq_dict = {}
    for residue in chain.get_residues():
        res_num = int(residue.id[1])
        aminoacid = seq1(residue.resname, custom_map=protein_letters_3to1)
        seq_dict[res_num] = aminoacid
    pdb_sequence = [seq_dict[x] for x in sorted(seq_dict)]
    return ''.join([x for x in pdb_sequence])

def match_pdb_residue_num_to_seq(chain, ref=None):
    """Match PDB residue numbering (as given in PDB file) to
    a reference sequence (can be pdb sequence) numbered by index.

    Reference sequence is 1-indexed (and is indexed as such in output).

    Args:
        chain: A Bio.PDB chain object.
        ref (str): A reference protein sequence. Defaults to protein sequence
            given in PDB file.
    Returns:
        dict: A dictionary mapping reference sequence index (key) to
            residue numbering as given in the PDB file (value).
    """
    if ref is None:
        ref = chain.sequence
    seq_dict = {}
    for residue in chain.chain.get_residues():
        res_num = int(residue.id[1])
        aminoacid = seq1(residue.resname, custom_map=protein_letters_3to1)
        seq_dict[res_num] = aminoacid
    pdb_sequence = [seq_dict[x] for x in sorted(seq_dict)]
    pdb_numbering = [x for x in sorted(seq_dict)]
    pdb_to_ref, _ref_to_pdb = blast_sequences(pdb_sequence, ref)
    output = {}
    for i, pdb_num in enumerate(pdb_numbering):
        if i+1 in pdb_to_ref:
            output[pdb_to_ref[i+1]] = pdb_num
    return output

def _count_residues(_chain, _data, residues, _ref):
    """Simple function to count the number of residues within a radius"""
    return len(residues)

def _tajimas_d(_chain, alignment, residues, ref):
    """"Calculate Tajimas D for selected residues within a PDB chain.
    input is Chain object, multiple sequence alignment object,
    list of surrounding residues, and a dictionary giving mapping
    of PDB residue number to codon positions.
    """
    #filter list of residues based on those that have mapped codons:
    residues = [x for x in residues if x in ref]
    #Get list of codons that correspond to selected residues
    codons = [ref[res] for res in residues]
    #Get alignment bp from selected codons
    sub_align = _construct_sub_align(alignment, codons, fasta=True)
    #Compute Tajima's D using selected codons.
    tajd = gentests.tajimas_d(sub_align)
    return tajd

def _default_mapping(_chain, data, residues, ref):
    """"Calculate an average of all data points over selected residues.
    """
    #filter list of residues based on those that are mapped to reference seq
    residues = [x for x in residues if x in ref]
    #Convert PDB residue numbering to reference numbering
    reference_residues = [ref[res] for res in residues]
    data_points = [data[res] for res in reference_residues]
    average = np.mean(data_points)
    return average

def _snp_mapping(_chain, data, residues, ref):
    """"Calculate the percentage of SNPs over selected residues.
    Data is a list of residues that contain SNPs.
    """
    #filter list of residues based on those that are mapped to reference seq
    residues = [x for x in residues if x in ref]
    #Convert PDB residue numbering to reference numbering
    reference_residues = [ref[res] for res in residues]
    #Find the intersection between the residues which contain SNPs and
    #the selected residues on the Structure
    snp_xor_res = set(data) & set(reference_residues)
    num_snps = len(snp_xor_res)
    try:
        perc_snps = num_snps / len(reference_residues) * 100
    #If no residues are mapped onto the reference sequence, return None.
    except ZeroDivisionError:
        return None
    #currently returns the proportion of SNPs within a radius. could
    #change to be the raw number of SNPs.
    return perc_snps

def _map_amino_acid_scale(chain, data, residues, _ref):
    """
    Compute average value for amino acid propensity scale.
    """
    #Get a list of all amino acids within window, converted to one letter code
    aminoacids = [seq1(chain[int(res)].resname, custom_map=protein_letters_3to1)
                  for res in residues]
    scales = {'kd': ProtParamData.kd, # Kyte & Doolittle index of hydrophobicity
              # Flexibility
              # Normalized flexibility parameters (B-values),
              # average (Vihinen et al., 1994)
              'Flex': ProtParamData.Flex,
              # Hydrophilicity
              # Hopp & Wood
              # Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).
              'hw': ProtParamData.hw,
              # Surface accessibility
              # 1 Emini Surface fractional probability
              'em': ProtParamData.em,
              # 2 Janin Interior to surface transfer energy scale
              'ja': ProtParamData.ja}
    if data in scales:
        scale = scales[data]
    else:
        scale = data
    #Compute mean of scale over all residues within window
    result = np.mean([scale[aa] for aa in aminoacids])
    return result
