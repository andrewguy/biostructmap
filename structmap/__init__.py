"""StructMap is a package for performing analyses on Protein Data Bank
(PDB) structure files. In particular, it allows for the spatial averaging
of sequence-aligned data over a given protein structure, and also allows the
incorporation of structural information into genetic tests of selection pressure
such as Tajima's D.
"""

from .structmap import Structure, SequenceAlignment

__version__ = '0.2.0'
__all__ = ["structmap", "seqtools", "pdbtools", "gentests", "map_functions"]
