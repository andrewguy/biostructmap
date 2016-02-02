import Bio.PDB
from structmap import utils, pdbtools, seqtools, gentests
from Bio import SeqIO, AlignIO
from Bio.PDB import DSSP

class Structure:
    def __init__(self, pdbfile):
        #create pdb parser, get structure.
        parser = Bio.PDB.PDBParser()
        pdbname = pdbfile
        #Get Bio.PDB structure
        self.structure = parser.get_structure(pdbname,pdbfile)
        #Get PDB sequences as a list of strings
        self.sequence = pdbtools.get_pdb_seq(pdbfile)
        self.chains = {}
        for model in self.structure:
            for chain in model:
                self.chains[chain.get_id()]=Chain(chain, pdbfile)

class Chain:
    def __init__(self, chain,pdbfile):
        self.chain = chain
        model = self.chain.get_parent.get_id()
        self.dssp = DSSP(model,pdb_file)

    def nearby(self,radius=15,atom='all'):
        '''Take a Bio.PDB chain object, and find all residues within a radius of
        a given residue. Return a dictionary containing nearby residues for each
        residue in the chain.
        Optional parameter is the atom with which to compute distance.
        By default this is 'all', which gets all non-heterologous atoms.
        Other potential options include 'CA', 'CB' etc. If an atom is not found
        within a residue object, then method reverts to using 'CA'.
        '''
        self.dist_map = pdbtools.nearby(self.chain,radius,atom)
        return self.dist_map

    def rsa(self):
        '''Use Bio.PDB to calculate relative solvent accessibility.
        Return a dictionary with RSA values for each residue.
        '''
        pass


class Sequence:
    def __init__(self, seqfile):
        self.seqrecord = SeqIO.read(open(seqfile), "fasta")

    def sequence(self):
        sequence = utils.to_string(self.seqrecord)
        return sequence


class SequenceAlignment:
    def __init__(self, alignfile, file_format='fasta'):
        self.alignment = AlignIO.read(alignfile,file_format)

    def tajimas_d(self,window=None,step=3):
        return gentests.tajimas_d(self.alignment,window,step)
