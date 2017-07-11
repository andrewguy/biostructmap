"""
Main PDB Structure Mapping (StructMap) module. This package implements
techniques to map data onto a protein PDB structure.

Rationale:

Whilst many protein-related data (such as genetic polymorphisms) can be easily
related to a reference protein sequence, it is slightly more involved to map
this data onto a protein structure. For example, when viewing polymorphisms
that have arisen as a result of immune selection pressure, it is instructive
to be able to view the location of these polymorphisms in the context of the
protein 3D structure. Furthermore, it may also be instructive to average this
data over a radius around a central residue. In our example case of protein
polymoprhisms, this allows the user to identify polymorphic hotspots on the
protein structure.

We can extend this concept further, and consider genetic tests for phenomena
which arise as a result of selection pressure at the level of protein structure.
For example, antibody-mediated immune selection pressure gives rise to balancing
selection at the level of population genetics. A number of methods exist to
identify regions of the genome under balancing selection, including Tajima's D.
Tajima's D is often calculated as a sliding window over a genome. However, when
the origins of balancing selection occur at the level of protein structure,
it can be useful to also consider protein spatial information when computing
Tajima's D using a sliding window. In essence, we could calculate Tajima's D
using a 3D sliding window over a protein structure. The StructMap package
automates this process.

Details:

The StructMap package makes extensive use of the Biopython Bio.PDB module for
PDB file parsing and integration with calculation of secondary structure
and relative solvent accessibility via the DSSP program. Calculation of Tajima's
D is performed using DendroPy, although several optimisations are performed
before passing a multiple sequence alignment to DendroPy in order to speed
up calculation of Tajima's D, such as removal of non-polymorphic sites and
memoization of results from previous windows when calculating a sliding window
value of Tajima's D.


"""
from __future__ import absolute_import, division, print_function

from copy import deepcopy
from os import path
from tempfile import NamedTemporaryFile
from Bio.PDB import DSSP, PDBIO, PDBParser, FastMMCIFParser
from Bio import AlignIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from . import pdbtools, gentests
from .pdbtools import match_pdb_residue_num_to_seq, SS_LOOKUP_DICT, mmcif_sequence_to_res_id
from .map_functions import (_tajimas_d, _default_mapping, _snp_mapping,
                            _map_amino_acid_scale)
from .seqtools import (blast_sequences, align_protein_to_dna,
                       _construct_sub_align)

class DataMap(dict):
    '''
    A class to hold a mapping of data to some PDB object. Extends the `dict`
    class, adding a few methods to allow correspondence between data and PDB
    chain object.

    Should be called as `DataMap(results, chain=chain, params=params)`,
    where `results` would be a dictionary of mapped data over a PDB chain for a
    particular set of parameters.
    '''
    def __init__(self, *args, **kw):
        """Initialise a DataMap object, which stores mapped to some PDB object
        as a dictionary, but also provides a link to original chain object
        and parameters used to generate data.

        Args:
            *args: Standard dict args
            **kwargs: Standard dict kwargs. Kwargs `chain` and `params` are
                reserved for reference to the PDB chain object and analysis
                parameters respectively, and are required. These are removed
                from the list of kwargs before passing to dict __init__ method.
        """
        self.chain = kw.pop('chain')
        self.params = kw.pop('params')
        super(DataMap, self).__init__(*args, **kw)


    def write_data_to_pdb_b_factor(self, default_no_value=0, outdir='',
                                   filename=None):
        """Write mapped data to PDB B-factor column, and save as a PDB file.

        This method allows for data visualisation over the PDB structure
        using Pymol. In particular, the `spectrum` command in Pymol can color
        residues according to B-factor values, and it is common practice to
        store data of interest in the B-factor column of PDB files.

        Args:
            default_no_value (int, optional): Default value to write to B-factor
                column if a non-numeric value exists for a particular data point
                (e.g. a None value). Default value is 0, although this could be
                altered depending the range of possible data points. For
                example, you might use a value such as 99 or -99 to indicate a
                missing data point, and then filter this out when visualising
                the data in Pymol.
            outdir (str, optional): Output directory to write PDB file to.
                Defaults to current directory.
            filename (str, optional): Output file name. If not provided,
                this defaults to a sensible description of the data being
                written, in the form of "pdbid_param1-value1_param2-value2.pdb".

        Returns:
            None
        """
        # Make copy of chain object. This may fail in earlier versions of
        # Biopython if chain contains disorded residues. We do this so that
        # we are not overwriting B-factor values in original PDB chain.
        _chain = deepcopy(self.chain)
        _structure = _chain.parent().parent()
        #Set all B-factor fields to zeros/default value
        for residue in _chain:
            _data = self.get(residue.get_id(), default_no_value)
            if _data is None:
                _data = default_no_value
            for atom in residue:
                if atom.is_disordered():
                    for altloc in atom.disordered_get_id_list():
                        atom.disordered_select(altloc)
                        atom.set_bfactor(float(_data))
                else:
                    atom.set_bfactor(float(_data))
        pdb_io = PDBIO()
        pdb_io.set_structure(_chain.chain)
        # Write chain to PDB file
        if filename is None:
            filename = _structure.pdbname + '_' + self._parameter_string() + '.pdb'
        pdb_io.save(path.join(outdir, filename))
        return None

    def write_to_atom(self, output, sep=','): #TODO write tests
        """Write score for each atom in a structure to a file, based on
        a data dictionary mapping output score to residue number.

        Each line of the output file contains an entry for a single atom.

        Args:
            output (str): Output file name/path.
            sep (str, optional): Seperator between residue and data.
                Defaults to `,`.
        Returns:
            None
        """
        #For each atom in the structure, write an output score based on the data
        #given, presuming (key,value) in a dictionary with key corresponding to
        #a residue number.
        with open(output, 'w') as f:
            for res in sorted(self):
                residue = self.chain[int(res)]
                for atom in residue:
                    data_pt = [str(x) for x in [atom.serial_number, self[res]]]
                    line = sep.join(data_pt) + '\n'
                    f.write(line)
        return

    def write_to_residue(self, output, sep=',', ref=None): #TODO write tests
        """Write score for each residue in a structure to a file, based on
        a data dictionary mapping output score to residue number.

        Each line of the output file contains an entry for a single residue.
        This method additionally allows for scores to be aligned to a reference
        sequence for easy comparison with other metrics or a conventional
        reference sequence.

        Args:
            output (str): Output file name/path.
            sep (str, optional): Seperator between residue and data.
                Defaults to `,`.
            ref (str, optional): A reference sequence with which to align
                data to.
        Returns:
            None
        """
        if ref is None:
            with open(output, 'w') as f:
                for res in sorted(self):
                    data_pt = [str(x) for x in [res, self[res]]]
                    line = sep.join(data_pt) + '\n'
                    f.write(line)
        else:
            seq_index_to_pdb_numb = match_pdb_residue_num_to_seq(self.chain, self.chain.sequence)
            pdbindex_to_ref, _ = blast_sequences(self.chain.sequence, ref)
            pdbnum_to_ref = {seq_index_to_pdb_numb[x]:pdbindex_to_ref[x] for x
                             in pdbindex_to_ref if x in seq_index_to_pdb_numb}
            with open(output, 'w') as f:
                for res in sorted(self):
                    if res not in pdbnum_to_ref:
                        output = ("Residue {res} in PDB file {pdb} was not"
                                  " matched to reference sequence provided"
                                  " for writing to output file").format(
                                      res=res,
                                      pdb=self.chain.parent().parent().pdbname)

                        print(output)
                        continue
                    data_pt = [str(x) for x in [pdbnum_to_ref[res], self[res]]]
                    line = sep.join(data_pt) + '\n'
                    f.write(line)
        return

    def _parameter_string(self):
        """Create descriptive string from parameter values"""
        param_string = '_'.join(["{k}-{v}".format(k=key, v=value)
                                 for key, value in sorted(self.params.items())])
        return param_string


class Structure(object):
    """A class to hold a PDB structure object.

    Attributes:
        structure: Underlying Bio.PDB.Structure.Structure object
        sequences (dict): A dictionary of all protein sequences within
            structure, accessed by chain id.
        models (dict): Dictionary of all models within structure, accessed by
            model id.
        pdbname (str): A descriptive name for the PDB file.
    """
    def __init__(self, pdbfile, pdbname='pdb_file', mmcif=False):
        """Initialise PDB Structure object.

        Args:
            pdbfile (str/file-object): A filename string or a file-like object
                that contains a PDB file.
            pdbname (str): A descriptive name for the PDB file. This is used
                as part of the default naming options when writing output files.
            mmcif (bool, optional): Set to true if reading a PDBx/mmCIF file,
                otherwise this defaults to reading a PDB file.
        """
        # Create PDB parser
        if mmcif:
            parser = FastMMCIFParser()
        else:
            parser = PDBParser()
        self._mmcif_dict = None
        self._mmcif = mmcif
        #Get Bio.PDB structure
        self._pdbfile = pdbfile
        self.structure = parser.get_structure(pdbname, self.pdb_file())
        #Get PDB sequences
        if mmcif:
            self.sequences = pdbtools.get_mmcif_canonical_seq(self.mmcif_dict())
        else:
            self.sequences = pdbtools.get_pdb_seq(self.pdb_file())
        self.models = {model.get_id():Model(self, model) for
                       model in self.structure}
        self.pdbname = pdbname


    def __iter__(self):
        """Iterate over all models within structure"""
        for key in sorted(self.models):
            yield self.models[key]

    def __getitem__(self, key):
        """Get a model from structure"""
        return self.models[key]

    def pdb_file(self):
        '''Return the PDB file object, which can either be a string, or a
        file-like object.

        Notes:
            If it is a file-like object, sets read point to start of the
            file or stream. This enables us to perform analysis on PDB files
            downloaded from RCSB and stored in memory (as a StringIO object).

        Returns:
            str/file object: Either a string representing a filename or a file
                object with read point set to the start of the file or stream.
        '''
        if not isinstance(self._pdbfile, str):
            self._pdbfile.seek(0)
        return self._pdbfile

    def mmcif_dict(self):
        '''Return the mmcif dictionary.

        Only applicable if using an mmCIF file.

        Returns:
            dict: A dictionary containing mmCIF data.
        '''
        if self._mmcif and self._mmcif_dict is None:
            self._mmcif_dict = MMCIF2Dict(self.pdb_file())
        elif not self._mmcif:
            raise TypeError("Not an mmCIF file!")
        return self._mmcif_dict


class Model(object):
    """A class to hold a PDB model object.

    Attributes:
        model: Underlying Bio.PDB.Model.Model object
        dssp (dict): DSSP results for the model, ONLY IF model is the first in
            PDB file. Otherwise set to an empty dict. This is because the DSSP
            program will only read the first model in a PDB file.
        chains (dict): Dictionary of all chains within model, accessed by chain
            id.
    """
    def __init__(self, structure, model):
        """Initialise A PDB Model object.

        Args:
            structure (structmap.Structure): Parent structure object.
            model (Bio.PDB.Model.Model): Bio.PDB Model object.
        """
        self._id = model.get_id()
        self.model = model
        self._parent = structure
        #DSSP only works on the first model in the PDB file
        if isinstance(structure.pdb_file(), str) and self._id == 0:
            try:
                self.dssp = DSSP(self.model, structure.pdb_file())
            except OSError:
                self.dssp = DSSP(self.model, structure.pdb_file(),
                                 dssp="mkdssp")
        elif self._id == 0:
            if structure._mmcif:
                suffix = '.cif'
            else:
                suffix = '.pdb'
            with NamedTemporaryFile(mode='w', suffix=suffix) as temp_pdb_file:
                temp_pdb_file.write(structure.pdb_file().read())
                temp_pdb_file.flush()
                try:
                    self.dssp = DSSP(self.model, temp_pdb_file.name)
                except OSError:
                    self.dssp = DSSP(self.model, temp_pdb_file.name,
                                     dssp="mkdssp")
        else:
            self.dssp = {}
        self.chains = {chain.get_id():Chain(self, chain) for
                       chain in self.model}

    def __iter__(self):
        """Iterate over all chains in model"""
        for key in sorted(self.chains):
            yield self.chains[key]

    def __getitem__(self, key):
        """Get selected chain"""
        return self.chains[key]

    def parent(self):
        """Get parent Structure object"""
        return self._parent

    def get_id(self):
        """Get model ID"""
        return self._id


class Chain(object):
    """A class to hold a PDB chain object.

    Attributes:
        chain: Underlying Bio.PDB.Chain.Chain object
        dssp (dict): DSSP results for the parent model, ONLY IF model is the
            first in PDB file. Otherwise set to an empty dict. This is because
            the DSSP program will only read the first model in a PDB file.
        sequence (str): Chain protein sequence. Extracted from PDB file header
            if present, otherwise constructed from PDB atom records for that
            chain.
    """
    def __init__(self, model, chain):
        """Initialise A PDB Chain object.

        Args:
            model (structmap.Model): Parent Model object.
            chain (Bio.PDB.Chain.Chain): Bio.PDB Chain object.
        """
        self._id = chain.get_id()
        self.chain = chain
        self._parent = model
        #Refer to model for DSSP data
        self.dssp = self._parent.dssp
        #Some PDB files do not contain sequences in the header, and
        #hence we need to parse the atom records for each chain
        try:
            self.sequence = model.parent().sequences[self.get_id()]
        except KeyError:
            if self._parent._parent._mmcif:
                raise KeyError("mmCIF file doesn't contain sequences " \
                              "for chain {chain_id}".format(chain_id=self._id))
            self.sequence = pdbtools.get_pdb_seq_from_atom(chain)
            self._parent.parent().sequences[self.get_id()] = self.sequence
        self._nearby = {}
        self._rsa = {}

    def __iter__(self):
        """Iterate over all residues in the chain"""
        for residue in self.chain:
            yield residue

    def __getitem__(self, key):
        """Get residue within chain"""
        return self.chain[key]

    def parent(self):
        """Get parent Model object"""
        return self._parent

    def get_id(self):
        """Get chain ID"""
        return self._id

    def nearby(self, radius=15, atom='all'):
        """Take a Bio.PDB chain object, and find all residues within a radius of
        a given residue.

        Args:
            radius (int/float): Radius within which to find nearby residues for
                each residue in the chain.
            atom (str): The atom with which to compute distances. By default
                this is 'all', which gets all non-heterologous atoms. Other
                potential options include 'CA', 'CB' etc. If an atom is not
                found within a residue object, then method reverts to using
                'CA'.

        Returns:
            dict: A dictionary containing a list of nearby residues for each
                residue in the chain.
        """
        param_key = (radius, atom)
        #Calculate distance matrix and store it for retrieval in future queries.
        if param_key not in self._nearby:
            dist_map = pdbtools.nearby(self.chain, radius, atom)
            self._nearby[param_key] = dist_map
        return self._nearby[param_key]

    def rel_solvent_access(self):
        """Use Bio.PDB DSSP tools to calculate relative solvent accessibility.

        Returns:
            dict: A dictionary with RSA values for each residue (key).
        """
        if not self._rsa:
            rsa = {}
            for residue in self.chain:
                key = (self.get_id(), residue.get_id())
                if key in self.dssp:
                    try:
                        rsa[key[1]] = float(self.dssp[key][3])
                    except ValueError:
                        rsa[key[1]] = None
            self._rsa = rsa
        return self._rsa

    def secondary_structure(self, numeric_ss_code=False):
        '''Use DSSP to calculate secondary structure elements.

        Notes:
            Secondary structure assignment notation follows that of the DSSP
            program.

            The DSSP codes for secondary structure used here are:
                =====     ====
                Code      Structure
                =====     ====
                 H        Alpha helix (4-12)
                 B        Isolated beta-bridge residue
                 E        Strand
                 G        3-10 helix
                 I        Pi helix
                 T        Turn
                 S        Bend
                 -       None
                =====     ====
            If kwarg `numeric_ss_code` == True, then secondary structure values
            will be returned as integers according to the following lookup
            table:
            {
            'H': 0,
            'B': 1,
            'E': 2,
            'G': 3,
            'I': 4,
            'T': 5,
            'S': 6,
            '-', 7
            }
            These are provided in pdbtools.SS_LOOKUP_DICT, as well as the
            reverse lookups (i.e. {0: 'H', ...}).

        Args:
            numeric_ss_code (bool): Return secondary structure elements as a
                numeric lookup code if True. Set to False by default.

        Returns:
            dict: A dictionary with secondary structure assignment (value)
                for each residue (key) within the chain.
        '''
        keys = [x for x in self.dssp.keys() if x[0] == self._id]
        ss_dict = {x[1]: self.dssp.property_dict[x][2] for x in keys}
        if numeric_ss_code:
            return {key:SS_LOOKUP_DICT[item] for key, item in ss_dict.items()}
        else:
            return ss_dict

    def map(self, data, method='default', ref=None, radius=15, selector='all',
            rsa_range=None, map_to_dna=False):
        """Perform a mapping of some parameter or function to a pdb structure,
        with the ability to apply the function over a '3D sliding window'.

        The residues within a radius of a central residue are passed to the
        mapping function, which computes an output value for the central
        residue. This is performed for each residue in the structure.
        For calculation of Tajima's D, reference sequence is a genomic sequence.

        Args:
            data (dict/object): Generally a dictionary containing data to
                to be mapped over structure. The exact form for this data
                depends on the function being used to map data - each mapping
                function is passed this data object.
            method (str): A string representing a method for mapping data.
                Internally, these strings are keys for a dictionary of
                functions, and these can be extended with custom user-provided
                functions if desired.
            ref (str): A reference protein sequence. For calculation of Tajima's
                D, this is a reference genomic sequence. Data to be mapped
                should be aligned/referenced according to this reference
                sequence.
            radius (int/float, optional): The radius (Angstrom) over which to
                select nearby residues for inclusion within each 3D window.
                This defaults to 15 Angstrom, which is the typical maximum
                dimension for an antibody epitope.
            selector (str, optional): A string indicating the atom with which
                to compute distances between residues. By default this is 'all',
                which gets all non-heterologous atoms. Other potential options
                include 'CA', 'CB' etc. If an atom is not found within a residue
                object, then the selection method reverts to using 'CA'.
            rsa_range (tuple, optional): A tuple giving (minimum, maximum)
                values of relative solvent accessibility with which to filter
                all residues on (ie. map method will ignore all residues
                outside this range). This is useful when wanting to examine only
                surface exposed residues.
            map_to_dna (bool, optional): Set True if the mapping method involves
                aligning to a DNA sequence. Defaults to False.

        Returns:
            structmap.DataMap: A dictionary-like object which contains mapped
                values for each residue (key). This object extends the standard
                dict type, adding methods to allow writing of data to PDB
                B-factor columns for easy viewing using Pymol or other similar
                programs.
        """
        #Note: This method attempts to deal with 3 different ways of identifying
        #residue position: i) Within a PDB file, residues are labelled with a
        #number, which should increment according to position in sequence, but
        #is not necessarily gapless, and may include negative numbering in
        #order to preserve alignment with similar sequences. ii) Residue
        #numbering according to the position of the residue within the
        #protein sequence extracted from the PDB file. That is, the first
        #residue in the sequence is numbering '1', and incrementing from there.
        #iii) Residue numbering according to a reference sequence provided, or
        #according to a dna sequence (which could include introns).

        methods = {"default":_default_mapping,
                   "tajimasd":_tajimas_d,
                   "snps": _snp_mapping,
                   "aa_scale": _map_amino_acid_scale}
        #Create a map of pdb sequence index (1-indexed) to pdb residue
        #numbering from file
        is_mmcif = self._parent._parent._mmcif
        if is_mmcif:
            mmcif_dict = self.parent().parent().mmcif_dict()
            _, _seq_index_to_pdb_numb = mmcif_sequence_to_res_id(mmcif_dict)
            seq_index_to_pdb_numb = {key: value[1] for key, value in
                                     _seq_index_to_pdb_numb.items() if
                                     value[0] == self.get_id()}
        else:
            seq_index_to_pdb_numb = match_pdb_residue_num_to_seq(self, self.sequence)
        #Generate a map of nearby residues for each residue in pdb file. numbering
        #is according to pdb residue numbering from file
        residue_map = self.nearby(radius=radius, atom=selector)
        results = {}

        # If method involves mapping to a genome, then we presume the data given
        # contains a multiple sequence alignment, and we use the first sequence
        # as reference. Note that it would be better for the user to explicitly
        # provide a reference genome in most cases.
        if ref is None and map_to_dna:
            ref = data[0]
        elif ref is None:
            ref = self.sequence
        #Generate mapping of pdb sequence index to dna sequence
        if map_to_dna:
            pdbindex_to_ref = align_protein_to_dna(self.sequence,
                                                   ''.join([x for x in ref]))
        #Generate mapping of pdb sequence index to reference sequence (also indexed by position)
        else:
            pdbindex_to_ref, _ = blast_sequences(self.sequence, ref)
        if method in methods:
            method = methods[method]
        #Finally, map pdb numbering by file to the reference sequence
        #(dna or protein) provided, as long as the residues exists within the PDB
        #structure (ie has coordinates)
        pdbnum_to_ref = {seq_index_to_pdb_numb[x]:pdbindex_to_ref[x] for x in
                         pdbindex_to_ref if x in seq_index_to_pdb_numb}
        results = {}
        #For each residue within the sequence, apply a function and return result.
        for residue in residue_map:
            if rsa_range:
                residues = self._filter_rsa(residue_map[residue], rsa_range)
                if residue not in residues:
                    results[residue] = None
                    continue
            else:
                residues = residue_map[residue]
            results[residue] = method(self, data, residues, pdbnum_to_ref)
        params = {'radius':radius, 'selector': selector}
        return DataMap(results, chain=self, params=params)

    def write_to_atom(self, data, output, sep=','): #TODO move to DataMap object
        """Write score for each atom in a structure to a file, based on
        a data dictionary mapping output score to residue number.

        Each line of the output file contains an entry for a single atom.

        Args:
            data (dict): A dictionary of output score for each residue. Can also
                be structmap.DataMap object, which extends the dict class.
            output (str): Output file name/path.
            sep (str, optional): Seperator between residue and data.
                Defaults to `,`.
        Returns:
            None
        """
        #For each atom in the structure, write an output score based on the data
        #given, presuming (key,value) in a dictionary with key corresponding to
        #a residue number.
        with open(output, 'w') as f:
            for res in data:
                residue = self.chain[int(res)]
                for atom in residue:
                    data_pt = [str(x) for x in [atom.serial_number, data[res]]]
                    line = sep.join(data_pt) + '\n'
                    f.write(line)
        return None

    def residue_to_atom_map(self):
        """Return a map of residue number (per numbering in PDB file) to atom
        number for this chain.

        Returns:
            dict: A dictionary mapping residue number (key) to a tuple (value)
                containing all atom serial numbers for that residue.
        """
        is_mmcif = self._parent._parent._mmcif
        if not is_mmcif:
            mapping = {residue.id[1]:tuple(atom.serial_number for atom in residue)
                       for residue in self.chain}
        else:
            raise TypeError("Can't map to atom serial with mmcif file!")
        return mapping

    def write_to_residue(self, data, output, sep=',', ref=None):
        """Write score for each residue in a structure to a file, based on
        a data dictionary mapping output score to residue number.

        Each line of the output file contains an entry for a single residue.
        This method additionally allows for scores to be aligned to a reference
        sequence for easy comparison with other metrics or a conventional
        reference sequence.

        Args:
            data (dict): A dictionary of output score for each residue. Can also
                be structmap.DataMap object, which extends the dict class.
            output (str): Output file name/path.
            sep (str, optional): Seperator between residue and data.
                Defaults to `,`.
            ref (str, optional): A reference sequence with which to align
                data to.
        Returns:
            None
        """
        if ref is None:
            with open(output, 'w') as f:
                for res in data:
                    data_pt = [str(x) for x in [res, data[res]]]
                    line = sep.join(data_pt) + '\n'
                    f.write(line)
        else:
            seq_index_to_pdb_numb = match_pdb_residue_num_to_seq(self, self.sequence)
            pdbindex_to_ref, _ = blast_sequences(self.sequence, ref)
            pdbnum_to_ref = {seq_index_to_pdb_numb[x]:pdbindex_to_ref[x] for x
                             in pdbindex_to_ref if x in seq_index_to_pdb_numb}
            with open(output, 'w') as f:
                for res in data:
                    if res not in pdbnum_to_ref:
                        output = ("Residue {res} in PDB file {pdb} was not"
                                  " matched to reference sequence provided"
                                  " for writing to output file").format(
                                      res=res,
                                      pdb=self.parent().parent().pdbname)

                        print(output)
                        continue
                    data_pt = [str(x) for x in [pdbnum_to_ref[res], data[res]]]
                    line = sep.join(data_pt) + '\n'
                    f.write(line)
        return None

    def _filter_rsa(self, residues, rsa_range):
        '''
        Function to remove residues with relative solvent accessibility
        values outside the requested range.

        If a residue solvent accessibility cannot be calculated by DSSP,
        then that residue is ignored.

        Args:
            residues (list): List of residues to filter.
            rsa_range (list/tuple): A tuple/list of the form (min, max),
                that gives the minimum and maximum relative solvent
                accessibility values with which to filter residues. Note that
                these values should be between 0 and 1.

        Returns:
            list: A list of residues (int) that have a relative solvent
                accessibility within the required range (inclusive).
        '''
        rsa = self.rel_solvent_access()
        filtered_residues = [x for x in residues if rsa.get(x, None) and
                             rsa_range[0] <= rsa[x] <= rsa_range[1]]
        return filtered_residues


class SequenceAlignment(object):
    """A class to hold a multiple sequence alignment object.

    This class is a basic wrapper for a Bio.Align.MultipleSequenceAlignment
    object, and allows calculation of Tajima's D over selected codons, as well
    as methods for slightly more efficient selection of sub-alignments based on
    selected codons, which may or may not be continuous.

    Attributes:
        alignment (Bio.Align.MultipleSequenceAlignment): Multiple Sequence
            Alignment object from Bio.Align.MultipleSequenceAlignment.
    """
    def __init__(self, alignfile, file_format='fasta'):
        """Initialises a SequenceAlignment object with a multiple sequence
        alignment file or file-like object.

        Args:
            alignfile (str/file handle): A multiple sequence alignment
                file path (str) or file-like object. This is the same as the
                Bio.AlignIO.read(...) inputs.
            file_format (str, optional): File format that the sequence
                alignment is stored as. This is the same as the arguments for
                Bio.AlignIO.read(...), and defaults to 'fasta'.
        """
        self.alignment = AlignIO.read(alignfile, file_format)
        self._alignment_position_dict = None
        self._isolate_ids = None

    def __getitem__(self, key):
        return self.alignment[key]

    def get_alignment_position_dict(self):
        """Returns a dictionary with single base pair alignments for each
        position in the multiple sequence alignment.

        This method utilises memoization to speed up multiple slicing operations
        on MultipleSeqAlignment objects, which can be expensive. This is
        patricularly useful when constructing a sub-alignment representing
        selected residues in a protein structure, as we would typically have
        to make many such sub-alignments when iterating over all residues within
        a structure.

        Returns:
            dict: A dictionary containing Bio.Align.MultipleSeqAlignment objects
                with a single base pair represented for each indexed position.
                Dictionary is of the form {int: MultipleSeqAlignment, ...}.
        """
        if self._alignment_position_dict is None:
            self._alignment_position_dict = {}
            for i in range(len(self.alignment[0])):
                self._alignment_position_dict[i] = self.alignment[:, i]

        return self._alignment_position_dict

    def get_isolate_ids(self):
        """Returns the isolate IDs for each sequence in the multiple sequence
        alignment.

        Returns:
            list: List of isolates/strains as strings.
        """
        if self._isolate_ids is None:
            self._isolate_ids = [seq.id for seq in self.alignment]
        return self._isolate_ids


    def tajimas_d(self, window=None, step=3, protein_ref=None, genome_ref=None,
                  output_protein_num=False):
        """Calculate Tajima's D on a SequenceAlignment object.

        If no window parameter is passed to the function, then the function
        calculates Tajima's D over the whole sequence alignment and
        returns a single numerical result.

        If a window size is given, then the function returns a dictionary
        of Tajima's D values over a sliding window.

        Note:

        Args:
            window (int, optional): The size of the sliding window over which
                Tajima's D is calculated. Default is None, in which case a
                single Tajima's D value is calculated for the multiple sequence
                alignment.
            step (int, optional): Step size for sliding window calculation.
                Default step size of 3 (ie. one codon).
            protein_ref (str, optional): A reference protein sequence which is
                used to optionally determine the protein coding region of the
                genome, and only calculate Tajima's D over these regions.
            genome_ref (str/Bio.Seq.Seq, optional): A genome reference sequence
                that is aligned to the multiple sequence alignment. Required if
                providing a protein reference sequence.
            output_protein_num (bool, optional): If set to True, will also
                return Tajima's D values for each position in the protein
                reference sequence.
        Returns:
            float/dict: If window parameter is None, returns a single value for
                Tajima's D. Otherwise a dict of genome window midpoint to
                calculated Tajima's D values is returned
            list (optional): Map of protein residue number to calculated
                Tajima's D values as tuple pairs of the form (residue, TajimasD)
                This is not returned as a dictionary as some residues may have
                multiple mapped Tajima's D value depending on step size.
        """
        # If given a protein reference sequence, align genome to reference
        # sequence, and only perform a Tajima's D analysis over the protein
        # coding region. Otherwise perform a Tajima's D test over the full
        # genomic sequence.
        if genome_ref is not None and protein_ref is not None:
            #Align genome to reference
            prot_to_genome = align_protein_to_dna(protein_ref, ''.join([x for x in genome_ref]))
            #Get sorted list of codons
            codons = [prot_to_genome[x] for x in sorted(prot_to_genome)]
            #Construct a sub-alignment
            alignment = _construct_sub_align(self, codons, fasta=True)
        elif genome_ref is None and protein_ref is None:
            alignment = self.alignment
        elif protein_ref is None:
            raise TypeError("Missing protein_ref assignment")
        elif genome_ref is None:
            raise TypeError("Missing genome_ref assignment")

        #Perform Tajima's D calculation
        try:
            tajd = gentests.tajimas_d(alignment, window, step)
        except TypeError:
            print("Error calculating Tajima's D. Please check inputs to " +
                  "Tajima's D function.")
            raise
        if window is not None and protein_ref is not None:
            #Unpack list of codons
            codon_list = [x for sublist in codons for x in sublist]
            #Reference to genome numbering
            tajd_ref_to_genome = {codon_list[int(x) - 1]: tajd[x] for x in tajd}
            #Dictionary comprehension to reverse and unpack map of protein
            #position to genome codons. Initial dictionary is in the form
            # {1:(1,2,3), 2:(4,5,6,),...}, hence the need for the multiple
            #layers within the dictionary comprehension (needed to unpack each
            #codon tuple).
            genome_to_prot = {i:x for x in prot_to_genome for i in prot_to_genome[x]}
            tajd_ref_to_protein = [(genome_to_prot[x], tajd_ref_to_genome[x])
                                   for x in tajd_ref_to_genome]
            if output_protein_num:
                return tajd_ref_to_genome, tajd_ref_to_protein
            else:
                return tajd_ref_to_genome
        else:
            return tajd
