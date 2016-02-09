"""Main Structmap module.
Does things such as mapping of Tajima's D to a protein structure.
This doc-string needs updating."""
from __future__ import absolute_import, division, print_function

import Bio.PDB
from Bio.PDB import DSSP
from Bio import SeqIO, AlignIO
from structmap import utils, pdbtools, gentests
from structmap.pdbtools import _tajimas_d, _default_mapping
from structmap.seqtools import map_to_sequence


class Structure(object):
    """A class to hold a PDB structure object."""
    def __init__(self, pdbfile):
        #create pdb parser, get structure.
        parser = Bio.PDB.PDBParser()
        pdbname = pdbfile
        #Get Bio.PDB structure
        self.structure = parser.get_structure(pdbname, pdbfile)
        #Get PDB sequences
        self.sequences = pdbtools.get_pdb_seq(pdbfile)
        self.models = {model.get_id():Model(self, model, pdbfile) for
                       model in self.structure}

    def __iter__(self):
        for key in sorted(self.models):
            yield self.models[key]

    def __getitem__(self, key):
        return self.models[key]

    def map(self, data, method='default', ref=None, radius=15, selector='all'):
        """Function which performs a mapping of some parameter or function to
        a pdb structure, with the ability to apply the function over a
        '3D sliding window'. The residues within a radius of a central
        residue are passed to the function, which computes an output value for
        the central residue. This is performed for each residue in the
        structure.
        """
        output = {}
        for model in self:
            for chain in model:
                unique_id = (model.get_id(), chain.get_id())
                output[unique_id] = chain.map(data, method=method, ref=ref,
                                              radius=radius, selector=selector)
        return output


class Model(object):
    """A class to hold a PDB model object. """
    def __init__(self, structure, model, pdbfile):
        self._id = model.get_id()
        self.model = model
        self._parent = structure
        self.chains = {chain.get_id():Chain(self, chain, pdbfile) for
                       chain in self.model}

    def __iter__(self):
        for key in sorted(self.chains):
            yield self.chains[key]

    def __getitem__(self, key):
        return self.chains[key]

    def parent(self):
        """Get parent Structure object"""
        return self._parent

    def get_id(self):
        """Get model ID"""
        return self._id


class Chain(object):
    """A class to hold a PDB chain object."""
    def __init__(self, model, chain, pdbfile):
        self._id = chain.get_id()
        self.chain = chain
        self._parent = model
        self.dssp = DSSP(model.model, pdbfile)
        self.sequence = model.parent().sequences[self.get_id()]

    def __iter__(self):
        for residue in self.chain:
            yield residue

    def __getitem__(self, key):
        return self.chain[key]

    def parent(self):
        """Get parent Model object"""
        return self._parent

    def get_id(self):
        """Get chain ID"""
        return self._id

    def nearby(self, radius=15, atom='all'):
        """Take a Bio.PDB chain object, and find all residues within a radius of
        a given residue. Return a dictionary containing nearby residues for each
        residue in the chain.
        Optional parameter is the atom with which to compute distance.
        By default this is 'all', which gets all non-heterologous atoms.
        Other potential options include 'CA', 'CB' etc. If an atom is not found
        within a residue object, then method reverts to using 'CA'.
        """
        dist_map = pdbtools.nearby(self.chain, radius, atom)
        return dist_map

    def rel_solvent_access(self):
        """Use Bio.PDB to calculate relative solvent accessibility.
        Return a dictionary with RSA values for each residue.
        """
        rsa = {}
        for residue in self.chain:
            key = (self.get_id(), residue.get_id())
            if key in self.dssp:
                rsa[key] = self.dssp[key][3]
        return rsa

    def map(self, data, method='default', ref=None, radius=15, selector='all'):
        """Perform a mapping of some parameter or function to a pdb structure,
        with the ability to apply the function over a '3D sliding window'.
        The residues within a radius of a central residue are passed to the
        function, which computes an output value for the central residue.
        This is performed for each residue in the structure.
        """
        methods = {"default":_default_mapping,
                   "tajimasd":_tajimas_d}
        residue_map = pdbtools.nearby(self.chain, radius=radius, selector=selector)
        results = {}
        if method == 'tajimasd' and ref is None:
            ref = data.translate(0)
        elif ref is None:
            ref = self.sequence
        pdb_to_ref, _ = map_to_sequence(self.sequence, ref)
        if method in methods:
            method = methods[method]
        for residue in residue_map:
            results[residue] = pdbtools.map_function(self, method, data,
                                                     residue_map[residue],
                                                     ref=pdb_to_ref)
        return results

    def write_to_atom(self, data, output):
        """Write score for each atom in a structure to a file, based on
        a data dictionary mapping output score to residue number.
        """
        #For each atom in the structure, write an output score based on the data
        #given, presuming (key,value) in a dictionary with key corresponding to
        #a residue number.
        with open(output, 'w') as f:
            for res in data:
                residue = self.chain[int(res)]
                for atom in residue:
                    data_pt = [str(x) for x in [atom.serial_number, data[res]]]
                    line = ','.join(data_pt) + '\n'
                    f.write(line)


class Sequence(object):
    """A class to hold a protein sequence"""
    def __init__(self, seqfile):
        self.seqrecord = SeqIO.read(open(seqfile), "fasta")

    def sequence(self):
        """Get protein sequence as a string"""
        sequence = utils.to_string(self.seqrecord)
        return sequence


class SequenceAlignment(object):
    """A class to hold a multiple sequence alignment object.
    Methods are:
        translate - translate to a protein sequence
        tajimas_d - calculate Tajima's D over a sliding window, or for the
        entire sequence.
    """
    def __init__(self, alignfile, file_format='fasta'):
        self.alignment = AlignIO.read(alignfile, file_format)

    def __getitem__(self, key):
        return self.alignment[key]

    def translate(self, index):
        """Translate to protein sequence. Translates until stop codon.
        Trims sequence to a multiple of 3.
        """
        length = len(self.alignment[index])
        overhang = length % 3
        translation = self.alignment[index,0:length-overhang].seq.translate(to_stop=True)
        return translation

    def tajimas_d(self, window=None, step=3):
        """Calculate Tajima's D on a SequenceAlignment object.

        If no window parameter is passed to the function, then the function
        calculates Tajima's D over the whole sequence alignment and
        returns a single numerical result.

        If a window size is given, then the function returns a dictionary
        of Tajima's D values

        :param window: The size of the sliding window over which Tajima's D is
        calculated
        :type window: int
        :param step: Step size for sliding window calculation
        :type step: int
        :returns: *key: window midpoint
                  *value: Tajima's D value for window
        """
        try:
            return gentests.tajimas_d(self.alignment, window, step)
        except TypeError:
            print("Error calculating Tajima's D. Please check inputs to " +
                  "Tajima's D function.")
            raise TypeError
