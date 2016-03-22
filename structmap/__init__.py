"""Structure mapper is a package for performing analyses on Protein Data Bank
(PDB) structure files.
"""

from .structmap import Structure, Sequence, SequenceAlignment

__version__ = '0.1.1'
__all__ = ["structmap", "seqtools", "pdbtools", "gentests", "utils"]
