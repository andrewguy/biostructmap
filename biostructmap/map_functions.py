'''A collection of methods that can be used to map data to a set of residues
within a PDB structure.

These methods all take the form:

def some_method(structure, data, residues, ref):
    ...
    return final_data

Required input parameters include:
    structure: a pdb structure object (Structure)
    data: a set of initial sequence-aligned data that is later filtered,
        and some function used to summarise filtered data.
    residues: a list of residues over which to filter input data. Note that
        filtering of initial data is not *required*, but is normally
        incorporated within the body of the defined function.
    ref: a dictionary mapping PDB residue numbers (key) to a reference sequence
        where the first residue in the sequence is indexed as 1, and
        incremented from there (value).

It is possible to write user-defined functions to pass to the Structure.map()
method following the above format ie. for the incorporation of novel amino
acid scales.

Note that departures from these input parameter conventions will be noted
in function docstrings, but otherwise can be assumed to follow the above
parameter convention.
'''
from __future__ import absolute_import, division

from Bio.SeqUtils import ProtParamData
from Bio.Data import IUPACData
import numpy as np
from .seqtools import _construct_sub_align_from_chains, _construct_protein_sub_align_from_chains
from . import gentests, protein_tests

IUPAC_3TO1_UPPER = {key.upper(): value for key, value in
                    IUPACData.protein_letters_3to1.items()}


def _count_residues(_structure, _data, residues, _ref):
    '''Simple function to count the number of residues within a radius.

    Returns:
        int: Number of residues within a radius'''
    return len(residues)

def _tajimas_d(_structure, alignments, residues, ref):
    '''Calculate Tajimas D for selected residues within a PDB chain.
    input is Chain object, multiple sequence alignment object,
    list of surrounding residues, and a dictionary giving mapping
    of PDB residue number to codon positions.

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to codon positions
            relative to the supplied multiple sequence alignment.
    Returns:
        float: Tajima's D value. Returns None if Tajima's D is undefined.
    '''
    tajd = _genetic_test_wrapper(_structure, alignments, residues, ref,
                                 gentests.tajimas_d)
    return tajd

def _wattersons_theta(_structure, alignments, residues, ref):
    '''Calculate wattersons theta for selected residues within a PDB chain.

    Input is Chain object, multiple sequence alignment object,
    list of surrounding residues, and a dictionary giving mapping
    of PDB residue number to codon positions.

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to codon positions
            relative to the supplied multiple sequence alignment.
    Returns:
        float: Wattersons theta. Returns None if wattersons theta
            is undefined.
    '''
    theta = _genetic_test_wrapper(_structure, alignments, residues, ref,
                                  gentests.wattersons_theta)
    return theta


def _shannon_entropy(_structure, alignments, residues, ref, table='Standard',
                     protein_letters=IUPACData.protein_letters, gap='-',
                     is_protein=False):
    '''Calculate Shannon entropy values for residues within a PDB chain.

    This should ideally be performed with radius = 0, so that the Shannon
    entropy of a single position is calculated. Otherwise, this method
    returns the mean Shannon entropy for all residues within the radius.

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to codon positions
            relative to the supplied multiple sequence alignment.
        table: A codon lookup table used by the Bio.Seq.translate() method.
            See BioPython docs for possible options.
        protein_letters (str, optional): String of all protein letters being
            used to define the amino acid alphabet. Defaults to standard 20
            amino acids. If another alphabet is used (if your sequence contains
            non-standard amino acid), then the maximum Shannon entropy values
            will change accordingly.
        gap (str): Character to denote a gap in the sequence alignment. Used
            when translating from DNA to protein sequence.
        is_protein (bool): Set to True if input sequence alignment is a protein sequence.
            Defaults to False, which means the input alignments should be DNA sequences.
    Returns:
        float: Shannon entropy, value between 0 (perfectly conserved) and
            4.322 (all residues are found at equal proportions at that position)
    '''
    if not is_protein:
        entropy = _genetic_test_wrapper(_structure, alignments, residues, ref,
                                        gentests.shannon_entropy, table=table,
                                        protein_letters=protein_letters, gap=gap)
    else:
        entropy = _protein_msa_wrapper(_structure, alignments, residues, ref,
                                       protein_tests.shannon_entropy_protein,
                                       protein_letters=protein_letters, gap=gap)
    return entropy

def _normalized_shannon_entropy(_structure, alignments, residues, ref,
                                table='Standard',
                                protein_letters=IUPACData.protein_letters, gap='-',
                                is_protein=False):
    '''Calculate normalized Shannon entropy values for residues within a PDB chain.

    This should ideally be performed with radius = 0, so that the normalized
    Shannon entropy of a single position is calculated. Otherwise, this method
    returns the mean Shannon entropy for all residues within the radius.
    The normalized Shannon entropy is the same as the Shannon entropy divided
    by the maximum possible entropy, and hence falls within the range [0, 1].

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to codon positions
            relative to the supplied multiple sequence alignment.
        table: A codon lookup table used by the Bio.Seq.translate() method.
            See BioPython docs for possible options.
        protein_letters (str, optional): String of all protein letters being
            used to define the amino acid alphabet. Defaults to standard 20
            amino acids. If another alphabet is used (if your sequence contains
            non-standard amino acid), then the maximum Shannon entropy values
            will change accordingly.
        gap (str): Character to denote a gap in the sequence alignment. Used
            when translating from DNA to protein sequence.
        is_protein (bool): Set to True if input sequence alignment is a protein sequence.
            Defaults to False, which means the input alignments should be DNA sequences.
    Returns:
        float: Normalized Shannon entropy, value between 0 (perfectly conserved)
           and 1 (all residues are found at equal proportions at that position)
    '''
    if not is_protein:
        entropy = _genetic_test_wrapper(_structure, alignments, residues, ref,
                                        gentests.shannon_entropy, table=table,
                                        protein_letters=protein_letters, gap=gap,
                                        normalized=True)
    else:
        entropy = _protein_msa_wrapper(_structure, alignments, residues, ref,
                                       protein_tests.shannon_entropy_protein,
                                       protein_letters=protein_letters, gap=gap,
                                       normalized=True)
    return entropy


def _nucleotide_diversity(_structure, alignments, residues, ref):
    '''Calculate nucleotide diversity for selected residues within a PDB chain.

    Input is Chain object, multiple sequence alignment object,
    list of surrounding residues, and a dictionary giving mapping
    of PDB residue number to codon positions.

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to codon positions
            relative to the supplied multiple sequence alignment.
    Returns:
        float: Nucleotide diversity. Returns None if nucleotide diversity
            is undefined.
    '''
    nucleotide_diversity = _genetic_test_wrapper(_structure, alignments,
                                                 residues, ref,
                                                 gentests.nucleotide_diversity)
    return nucleotide_diversity

def _dnds(_structure, alignments, residues, ref, treefile):
    '''Calculate dn/ds for selected residues within a PDB chain.

    Input is Chain object, multiple sequence alignment object,
    list of surrounding residues, a dictionary giving mapping
    of PDB residue number to codon positions, and a treefile to pass to PAML.

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to codon positions
            relative to the supplied multiple sequence alignment.
        treefile (str/file-like): A treefile to pass to PAML for calculation of
            dN/dS.
    Returns:
        float: dN/dS value.
    '''
    #Not currently implemented properly
    raise NotImplementedError
    chains = alignments.keys()
    #filter list of residues based on those that have mapped codons:
    ref_residues = [ref[x] for x in residues if x in ref]  # Returns [('A', (15,16,17)), ...]
    # Remove duplicate items (i.e. same data point, different chains) from
    # list of residues and set chain identifier to match alignment keys.
    codons = set([(chain, x[1]) for chain in chains
                  for x in ref_residues if x[0] in chain])
    #Get alignment bp from selected codons
    sub_align = _construct_sub_align_from_chains(alignments, codons, fasta=False)
    #Run external tool over sub alignment.
    with tempfile.NamedTemporaryFile(mode='w') as seq_file:
        seq_file.write(sub_align)
        seq_file.flush()
        process = subprocess.run(["/opt/bin/possum", "-f", "dnafasta", "-q", "-v", seq_file.name],
                                 stdout=subprocess.PIPE)
    try:
        output = float(process.stdout.decode().strip().split('\t')[-1])
    except ValueError:
        output = None
    return output

def _default_mapping(_structure, data, residues, ref, ignore_duplicates=True,
                     method=np.mean):
    '''Apply a data aggregation function over all data points over selected residues.

    Args:
        ignore_duplicates (bool): Ignore duplicate data points (i.e. from
            identical chains in a multi-chain structure).
        method (function): A data aggregation function. Should take a list of
            numeric values as input, and return a single numeric value. Default
            function is the arithmetic mean.

    Returns:
        float: Aggregated value for all data points within a radius.
    '''
    chains = data.keys()
    #filter list of residues based on those that have mapped reference
    ref_residues = [ref[x] for x in residues if x in ref]
    residues = [(chain, x[1]) for chain in chains
                for x in ref_residues if x[0] in chain]
    if ignore_duplicates:
        residues = set(residues)
    data_points = [data[res[0]][res[1]] for res in residues]
    if data_points:
        result = method(data_points)
    else:
        result = None
    return result


def _snp_mapping(_structure, data, residues, ref, ignore_duplicates=True,
                 output_count=False):
    '''Calculate the percentage of SNPs over selected residues.
    Data is a list of residues that contain SNPs.

    Args:
        data: List of residues that are polymorphic, aligned to reference seq.
        ignore_duplicates (bool): Ignore duplicate data points (i.e. on two
            separate chains) if True.
        output_percent (bool): Output the total number of SNPs within each
            radius if True. Default behaviour is to return a percentage of SNPs.

    Returns:
        float: Proportion of residues within a radius that are polymorphic.
    '''
    chains = data.keys()
    #filter list of residues based on those that have mapped reference
    ref_residues = [ref[x] for x in residues if x in ref]
    residues = [(chain, x[1]) for chain in chains
                for x in ref_residues if x[0] in chain]
    if ignore_duplicates:
        residues = set(residues)

    data = [(chain, snp) for chain, snps in data.items() for snp in snps]
    #Find the intersection between the residues which contain SNPs and
    #the selected residues on the Structure
    snp_xor_res = [residue for residue in residues if residue in data]
    num_snps = len(snp_xor_res)
    if output_count:
        output = num_snps
    else:
        try:
            output = num_snps / len(residues) * 100
        #If no residues are mapped onto the reference sequence, return None.
        except ZeroDivisionError:
            output = None
    return output


def _map_amino_acid_scale(structure, data, residues, _ref):
    '''
    Compute average value for amino acid propensity scale.

    Args:
        data (str): A string representing an amino acid propensity scale.
            Options are:
            'kd' -  Kyte & Doolittle index of hydrophobicity
            'Flex' - Flexibility, Normalized flexibility parameters (B-values)
            'hw' - Hydrophilicity, Hopp & Wood
            'em' - Surface accessibility, Emini Surface fractional probability
            'ja' - Surface accessibility, Janin Interior to surface transfer
                   energy scale

    Returns:
        float: Average propensity scale score over residues within a radius.
    '''
    first_model = sorted(structure.models)[0]
    #Get a list of all amino acids within window, converted to one letter code
    aminoacids = [IUPAC_3TO1_UPPER.get(
        structure[first_model][res[0]][res[1]].resname, 'X')
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

def _genetic_test_wrapper(_structure, alignments, residues, ref, genetic_test,
                          **kwargs):
    '''Helper function to generate a multiple sequence alignment from selected
    codons and pass to given function.

    Input is Chain object, multiple sequence alignment object,
    list of surrounding residues, a dictionary giving mapping
    of PDB residue number to codon positions, and function to pass subset of
    sequence alignment to. Used for genetic tests such as Tajima's D and
    Watterson's theta etc.

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to codon positions
            relative to the supplied multiple sequence alignment.
        genetic_test: A function that takes a multiple sequence alignment as
            input and returns a numeric value.
    Returns:
        float: Returned value from `genetic_test` function.
    '''
    chains = alignments.keys()
    #filter list of residues based on those that have mapped codons:
    ref_residues = [ref[x] for x in residues if x in ref]  # Returns [('A', (15,16,17)), ...]
    # Remove duplicate items (i.e. same data point, different chains) from
    # list of residues and set chain identifier to match alignment keys.
    codons = set([(chain, x[1]) for chain in chains
                  for x in ref_residues if x[0] in chain])
    #Get alignment bp from selected codons
    sub_align = _construct_sub_align_from_chains(alignments, codons, fasta=True)
    #Compute Tajima's D using selected codons.
    score = genetic_test(sub_align, **kwargs)
    return score

def _protein_msa_wrapper(_structure, alignments, residues, ref, msa_function,
                          **kwargs):
    '''Helper function to generate a multiple sequence alignment from selected
    residues and pass to given function.

    Input is Chain object, multiple sequence alignment object,
    list of surrounding residues, a dictionary giving mapping
    of PDB residue number to MSA residue positions, and function to pass subset of
    sequence alignment to. Used for calculations such as Shannon Entropy.

    Args:
        alignments (dict): A dictionary of multiple sequence alignments
            for each unique chain in the protein structure. Dictionary keys
            should be chain IDs.
        ref: A dictionary mapping PDB residue number to protein residue number
            relative to the supplied multiple sequence alignment.
        msa_function: A function that takes a multiple sequence alignment as
            input and returns a numeric value.
    Returns:
        float: Returned value from `genetic_test` function.
    '''
    chains = alignments.keys()
    #filter list of residues based on those that have mapped codons:
    ref_residues = [ref[x] for x in residues if x in ref]  # Returns [('A', (15,16,17)), ...]
    # Remove duplicate items (i.e. same data point, different chains) from
    # list of residues and set chain identifier to match alignment keys.
    residues = set([(chain, x[1]) for chain in chains
                  for x in ref_residues if x[0] in chain])
    #Get alignment bp from selected codons
    sub_align = _construct_protein_sub_align_from_chains(alignments, residues, fasta=True)
    #Compute something using selected alignment.
    score = msa_function(sub_align, **kwargs)
    return score
