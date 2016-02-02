from Bio.SeqIO import PdbIO
from scipy.spatial import distance
import numpy as np

def _euclidean_distance_matrix(chain, atom='all'):
    '''Compute the Euclidean distance matrix for all atoms in a pdb chain.
    Return the euclidean distance matrix and a reference list of all atoms
    in the model (positionally matched to the euclidean matrix).
    Optional parameter is the atom in each residue with which to compute a
    distance matrix. The default setting is 'all', which gets all
    non-heterologous atoms. Other potential options include 'CA', 'CB' etc.
    If an atom is not found within a residue object, then method reverts to
    using 'CA'.
    '''
    reference = []
    coords = []
    #Filter on non-HET atoms
    for residue in (x for x in chain if x.get_id()[0]==' '):
        #If selecting based on all atoms within residue
        if atom=='all':
            for atom in residue:
                coords.append(atom.get_coord())
                reference.append(atom.get_full_id()[3][1])
        #If measuring distance on particular atoms
        else:
            if atom in residue:
                x = atom
            #Revert to carbon alpha if atom is not found
            else:
                x = 'CA'
            coords.append(residue[x].get_coord())
            reference.append(residue[x].get_full_id()[3][1])
    #Convert to a np array, and compute Euclidean distance.
    coord_array = np.array(coords)
    euclid_mat = distance.pdist(coord_array,'euclidean')
    #Convert to squareform matrix
    euclid_mat = distance.squareform(euclid_mat)
    ref_array = np.array(reference)
    return euclid_mat, ref_array

def nearby(chain,radius=15,atom='all'):
    '''
    Take a Bio.PDB chain object, and find all residues within a radius of a
    given residue. Return a dictionary containing nearby residues for each
    residue in the chain.
    Optional parameter is the atom with which to compute distance. By default
    this is 'all', which gets all non-heterologous atoms. Other potential
    options include 'CA', 'CB' etc. If an atom is not found within a residue
    object, then method reverts to using 'CA'.
    '''
    #Setup variables
    ref_dict = {}
    euclidean_distance, reference = _euclidean_distance_matrix(chain,atom)
    within_radius = euclidean_distance<radius
    near_map = within_radius * reference
    #Iterate over all atoms in Euclidean distance matrix.
    for i in range(len(near_map)):
        if near_map[i][i] == 0:
            continue
        if near_map[i][i] not in ref_dict:
            ref_dict[near_map[i][i]] = near_map[i][np.nonzero(near_map[i])]
        else:
            ref_dict[near_map[i][i]] = np.append(ref_dict[near_map[i][i]],
                                        near_map[i][np.nonzero(near_map[i])])
    #Remove non-unique values.
    for key in ref_dict:
        ref_dict[key] = np.unique(ref_dict[key])
    return ref_dict


def _get_nearby(residue, chain, radius='10', atom='CA'):
    '''
    DEPRECATED
    Takes Bio.PDB residue and chain objects as input,and computes the residues
    which are within a certain radius of the central residue given. Output is a
    list containing the position of all nearby residues, including the central
    residue. Optional arguments include the radius within which to measure, and
    the particular atom from which to measure (carbon-alpha by default).
    '''
    nearby_residues = []
    for x in (x for x in chain if x.get_id()[0]==' '):
        if x[atom] - residue[atom] <= radius:
            nearby_residues.append(x.get_id()[1])
    return nearby_residues


def _map_surrounding_residues(chain, radius='10', atom='CA'):
    '''
    DEPRECATED
    Takes Bio.PDB chain and generates a map of surrounding residues for each
    residue in the chain. Output is a dictionary with each residue as a key,
    and a list of surrounding residues as a value.
    '''
    surr_res_map = {}
    for residue in (x for x in chain if x.get_id()[0]==' '):
        nearby_residues = get_nearby(residue,chain,radius,atom)
        surr_res_map[residue.get_id()[1]] = nearby_residues
    return surr_res_map


def get_pdb_seq(filename):
    '''
    Get a protein sequence from a PDB file.
    Will return multiple sequences if PDB file contains several chains.
    '''
    #Open PDB file and get sequence data
    with open(filename, 'r') as f:
        sequence_iter = [s for s in PdbIO.PdbSeqresIterator(f)]
        #A bit of manipulation to get Seq object into a string.
        sequences = [''.join([x for x in s]) for s in sequence_iter]
    return sequences
