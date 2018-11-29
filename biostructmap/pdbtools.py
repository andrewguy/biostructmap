"""A collection of tools for analysing a pdb file.

Helper module for the biostructmap package.
"""
from __future__ import absolute_import, division, print_function

from collections import defaultdict
from Bio.SeqIO import PdbIO
from Bio.SeqUtils import seq1
from Bio.Data.SCOPData import protein_letters_3to1
from Bio.PDB.Polypeptide import PPBuilder
import numpy as np
from .seqtools import align_protein_sequences
try:
    from scipy.spatial import distance, cKDTree
    SCIPY_PRESENT = True
except ImportError:
    SCIPY_PRESENT = False

SS_LOOKUP_DICT = {
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

def _euclidean_distance_matrix(model, selector='all'):
    """Compute the Euclidean distance matrix for all atoms in a pdb model.

    Args:
        model (Model): Bio.PDB Model object.
        selector (str): The atom in each residue with which to compute
            distances. The default setting is 'all', which gets all
            non-heterologous atoms. Other potential options include 'CA', 'CB'
            etc. If an atom is not found within a residue object, then method
            reverts to using 'CA'.
    Returns:
        np.array: A euclidean distance matrix.
        np.array: A reference list of all atoms in the model (positionally
            matched to the euclidean matrix).
    """
    reference = []
    coords = []
    # Get all non-HET residues from all chains
    residues = [res for chain in model for res in chain if
                res.get_id()[0] == ' ']
    #Filter on non-HET atoms
    for residue in residues:
        #If selecting based on all atoms within residue
        if selector == 'all':
            for atom in residue:
                coords.append(atom.get_coord())
                reference.append(atom.get_full_id()[2:4])
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
            reference.append(residue[select_atom].get_full_id()[2:4])
    #Convert to a np array, and compute Euclidean distance.
    coord_array = np.array(coords)
    euclid_mat = _pairwise_euclidean_distance(coord_array)
    ref_array = reference
    return euclid_mat, ref_array

def _pairwise_euclidean_distance(coord_array):
    '''Compute the pairwise euclidean distance matrix for a numpy array'''
    if SCIPY_PRESENT:
        euclid_mat = distance.pdist(coord_array, 'euclidean')
        #Convert to squareform matrix
        euclid_mat = distance.squareform(euclid_mat)
    else:
        euclid_mat = np.sqrt(((coord_array[:, :, None] -
                               coord_array[:, :, None].T) ** 2).sum(1))
    return euclid_mat

def _get_nearby_matrix(model, selector, radius):
    """Get a matrix of all nearby atoms in a pdb model.

    Args:
        model (Model): Bio.PDB Model object.
        selector (str): The atom in each residue with which to compute
            distances. The default setting is 'all', which gets all
            non-heterologous atoms. Other potential options include 'CA', 'CB'
            etc. If an atom is not found within a residue object, then method
            reverts to using 'CA'.
        radius (float): The radius within which to extract nearby residues/atoms.
    Returns:
        list: A nearby matrix (list of lists).
        np.array: A reference list of all atoms in the model (positionally
            matched to the nearby matrix).
    """
    reference = []
    coords = []
    # Get all non-HET residues from all chains
    residues = [res for chain in model for res in chain if
                res.get_id()[0] == ' ']
    #Filter on non-HET atoms
    for residue in residues:
        #If selecting based on all atoms within residue
        if selector == 'all':
            for atom in residue:
                coords.append(atom.get_coord())
                reference.append(atom.get_full_id()[2:4])
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
            reference.append(residue[select_atom].get_full_id()[2:4])
    #Convert to a np array, and use a KDTree to identify points within a certain distance.
    coord_array = np.array(coords)
    point_tree = cKDTree(coord_array)
    ball_tree = point_tree.query_ball_tree(point_tree, radius)
    ref_array = reference
    return ball_tree, ref_array

def nearby(model, radius=15, selector='all'):
    """
    Takes a Bio.PDB model object, and find all residues within a radius of a
    given residue.

    Args:
        model (Model): Bio.PDB Model object.
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
    #TODO: Cleanup case when scipy not present. Could just enforce scipy usage.
    ref_dict = {}
    #if SCIPY_PRESENT:
    if SCIPY_PRESENT:
        near_map, ref = _get_nearby_matrix(model, selector, radius)
        _ref_dict = defaultdict(set)
        for i, x in enumerate(near_map):
            _ref_dict[ref[i]].update({ref[y] for y in x})
    else:
        euclidean_distance, ref = _euclidean_distance_matrix(model, selector)
        within_radius = euclidean_distance <= radius
        del euclidean_distance
        # 1-indexed as 0 means not within range.
        near_map = within_radius * np.arange(1, len(ref)+1)
        #Iterate over all atoms in Euclidean distance matrix.
        for i, atom in enumerate(near_map):
            if atom[i] not in ref_dict:
                ref_dict[atom[i]] = atom[np.nonzero(atom)]
            else:
                ref_dict[atom[i]] = np.append(ref_dict[atom[i]],
                                              atom[np.nonzero(atom)])
        _ref_dict = {}
        del near_map
        # Go from numerical index to residue id
        for key, value in ref_dict.items():
            _ref_dict[ref[key-1]] = {ref[x-1] for x in value} | _ref_dict.get(ref[key-1], set())
    return _ref_dict


def mmcif_sequence_to_res_id(mmcif_dict):
    """Create a lookup from mmcif sequence id to a pdb residue ID and vice versa.

    This allows mapping between reference PDB sequences and BioPython
    residue ids.

    Args:
        mmcif_dict (dict): An mmcif dictionary from a Bio.PDB.Structure object.

    Returns:
        dict: Dictionary in the form {full Bio.PDB residue id: (chain,
              mmcif sequence id)}.
        dict: Dictionary in the form {(chain, mmcif sequence id):
              full Bio.PDB residue id}
    """
    _seq_id_list = mmcif_dict['_atom_site.label_seq_id']
    seq_id_list = []
    for seq_id in _seq_id_list:
        try:
            seq_id_list.append(int(seq_id))
        except ValueError:
            seq_id_list.append(None)
    # Parse dictionary to extract sequences from mmCIF file
    auth_chain_id_list = mmcif_dict['_atom_site.auth_asym_id']
    icode_list = mmcif_dict['_atom_site.pdbx_PDB_ins_code']
    # Create list of Bio.PDB full ids for each _atom_site
    hetero = mmcif_dict['_atom_site.group_PDB']
    residue_id_list = mmcif_dict['_atom_site.label_comp_id']
    # Use author provided numbering if given.
    if "_atom_site.auth_seq_id" in mmcif_dict:
        auth_seq_id_list = [int(x) for x in mmcif_dict["_atom_site.auth_seq_id"]]
    else:
        auth_seq_id_list = [int(x) for x in mmcif_dict["_atom_site.label_seq_id"]]
    # Get HET flags as used by Bio.PDB
    het_flag = [' ' if het == 'ATOM' else 'W' if res_type in ['HOH', 'WAT']
                else 'H_' + res_type for het, res_type in
                zip(hetero, residue_id_list)]
    # Get insertion code flags as used by Bio.PDB
    icode_flag = [' ' if icode == '?' else icode for icode in icode_list]
    # Residue ID as used by Bio.PDB
    bio_residue_ids = [*zip(het_flag, auth_seq_id_list, icode_flag)]
    full_ids = list(zip(auth_chain_id_list, bio_residue_ids))
    chain_and_seq_id_list = list(zip(auth_chain_id_list, seq_id_list))
    full_id_to_poly_seq_index = {key: value for key, value in
                                 zip(full_ids, chain_and_seq_id_list)}
    poly_seq_index_to_full_id = {value: key for key, value in
                                 full_id_to_poly_seq_index.items() if value}
    return full_id_to_poly_seq_index, poly_seq_index_to_full_id


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

def get_mmcif_canonical_seq(mmcif_dict):
    """
    Get structure sequences from an mmCIF file.

    Args:
        filename (str/filehandle): An mmCIF filename or file-like object.
    Returns:
        dict: Protein sequences (str) accesssed by chain id.
    """
    # TODO Should we use _pdbx_poly_seq_scheme here?
    # Parse dictionary to extract sequences from mmCIF file
    try:
        entity_seqs = mmcif_dict['_entity_poly.pdbx_seq_one_letter_code_can']
    except KeyError:
        entity_seqs = mmcif_dict['_entity_poly.pdbx_seq_one_letter_code']
    if isinstance(entity_seqs, list):
        chain_ids = [ids.split(',') for ids in
                     mmcif_dict['_entity_poly.pdbx_strand_id']]
        # Create dictionary of chain id (key) and sequences (value)
        sequences = dict((x, sublist[1].replace('\n', '')) for sublist in
                         zip(chain_ids, entity_seqs) for x in sublist[0])
    else:
        chain_ids = mmcif_dict['_entity_poly.pdbx_strand_id'].split(',')
        sequences = {chain_id: entity_seqs.replace('\n', '') for chain_id in chain_ids}
    return sequences

def get_mmcif_seqs(mmcif_dict):
    """
    Get structure sequences from an mmCIF file.

    Args:
        filename (str/filehandle): An mmCIF filename or file-like object.
    Returns:
        dict: Protein sequences (str) accesssed by chain id.
    """
    # TODO Should we use _pdbx_poly_seq_scheme here?
    # Parse dictionary to extract sequences from mmCIF file
    entity_ids = mmcif_dict['_entity_poly.entity_id']
    chain_ids = [ids.split(',') for ids in mmcif_dict['_entity_poly.pdbx_strand_id']]
    try:
        entity_seqs = mmcif_dict['_entity_poly.pdbx_seq_one_letter_code_can']
    except KeyError:
        entity_seqs = mmcif_dict['_entity_poly.pdbx_seq_one_letter_code']
    if isinstance(entity_seqs, list):
        # Create dictionary of chain id (key) and sequences (value)
        sequences = dict((x, sublist[1].replace('\n', '')) for sublist in
                         zip(chain_ids, entity_seqs) for x in sublist[0])
    else:
        sequences = {chain_ids[0]: entity_seqs.replace('\n', '')}
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
    # TODO Deprecate and revert to using polypeptide builder.
    seq_dict = {}
    for residue in chain.get_residues():
        res_num = int(residue.id[1])
        aminoacid = seq1(residue.resname, custom_map=protein_letters_3to1)
        seq_dict[res_num] = aminoacid
    # TODO Fix this, as not all residues are sequentially numbered.
    pdb_sequence = [seq_dict[x] for x in sorted(seq_dict)]
    return ''.join([x for x in pdb_sequence])

def match_pdb_residue_num_to_seq(model, ref=None):
    """Match PDB residue numbering (as given in PDB file) to
    a reference sequence (can be pdb sequence) numbered by index.

    Reference sequence is 1-indexed (and is indexed as such in output).

    Args:
        model: A biostructmap Model object.
        ref (dict): A dictionary containing reference protein sequences for each
            chain in the protein structure. Defaults to the protein sequences
            given in PDB file.
    Returns:
        dict: A dictionary mapping reference sequence index (key) to
            residue numbering as given in the PDB file (value). For example,
            we might have a key of ('A', 17) for the 17th residue in the
            reference sequence for chain 'A', with a value of
            ('A', (' ', 273, ' ')) that represents the Bio.PDB identifier for
            the corresponding residue.
    """
    ppb = PPBuilder()
    polypeptides = ppb.build_peptides(model.parent().structure)
    if ref is None:
        ref = model.parent().sequences
    output = {}
    for peptide in polypeptides:
        peptide_sequence = peptide.get_sequence()
        # Presume that peptide belongs to a single chain
        chain_id = peptide[0].get_full_id()[2]
        _, ref_to_pdb = align_protein_sequences(peptide_sequence, ref[chain_id])

        for ref_pos, pdb_pos in ref_to_pdb.items():
            output[(chain_id, ref_pos)] = peptide[pdb_pos - 1].get_full_id()[2:4]
    return output
