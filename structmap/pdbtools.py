"""Helper module for structmap package. Contains a collection of tools
for analysing a pdb file.
"""
from __future__ import absolute_import, division, print_function

from Bio.SeqIO import PdbIO
from scipy.spatial import distance
import numpy as np
from structmap.seqtools import (prot_to_dna_position,
                                _construct_sub_align)
from structmap import gentests

def _euclidean_distance_matrix(chain, selector='all'):
    """Compute the Euclidean distance matrix for all atoms in a pdb chain.
    Return the euclidean distance matrix and a reference list of all atoms
    in the model (positionally matched to the euclidean matrix).
    Optional parameter is the atom in each residue with which to compute a
    distance matrix. The default setting is 'all', which gets all
    non-heterologous atoms. Other potential options include 'CA', 'CB' etc.
    If an atom is not found within a residue object, then method reverts to
    using 'CA'.
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
            else:
                select_atom = 'CA'
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
    Take a Bio.PDB chain object, and find all residues within a radius of a
    given residue. Return a dictionary containing nearby residues for each
    residue in the chain.
    Optional parameter is the atom with which to compute distance. By default
    this is 'all', which gets all non-heterologous atoms. Other potential
    options include 'CA', 'CB' etc. If an atom is not found within a residue
    object, then method reverts to using 'CA'.
    """
    #Setup variables
    ref_dict = {}
    euclidean_distance, reference = _euclidean_distance_matrix(chain, selector)
    within_radius = euclidean_distance < radius
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
    """
    #Open PDB file and get sequence data
    with open(filename, 'r') as f:
        seq = [s for s in PdbIO.PdbSeqresIterator(f)]
        #A bit of manipulation to get Seq object into a dictionary
        #Key is chain ID, and value is sequence as a string.
        sequences = {s.id.split(":")[1]:''.join([x for x in s]) for s in seq}
    return sequences

def map_function(chain, method, data, residues, ref=None):
    """Map a function onto PDB residues, return an output value"""
    output = method(chain, data, residues, ref)
    return output

def _count_residues(_chain, _data, residues, _ref):
    """Simple function to count the number of residues within a radius"""
    return len(residues)

def _tajimas_d(_chain, alignment, residues, ref):
    """"Calculate Tajimas D for selected residues within a PDB chain.
    input is Chain object, multiple sequence alignment object,
    list of surrounding residues, and a dictionary giving mapping
    of PDB residue number to reference sequence residue number.
    """
    #Convert PDB residue numbering to reference numbering
    reference_residues = [ref[res] for res in residues]
    #Get lookup dictionary for bp position from residue numbers.
    codon_map = prot_to_dna_position(range(len(alignment[0])),
                                     range(len(alignment[0])//3))
    #Get list of codons that correspond to selected residues
    codons = [codon_map[res] for res in reference_residues]
    #Get alignment bp from selected codons
    sub_align = _construct_sub_align(alignment, codons)
    #Compute Tajima's D using selected codons.
    tajd = gentests.tajimas_d(sub_align)
    return tajd

def _default_mapping(_chain, data, residues, ref):
    """"Calculate an average of all data points over selected residues.
    """
    #Convert PDB residue numbering to reference numbering
    reference_residues = [ref[res] for res in residues]
    data_points = [data[res] for res in reference_residues]
    average = np.mean(data_points)
    return average
