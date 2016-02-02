from structmap import structmap

pdbfile = 'testpdb'
alignfile = 'testalign'

structure = structmap.Structure(pdbfile)
alignment = structmap.SequenceAlignment(alignfile)

tajimasd = alignment.tajimas_d()
taj_d_window = alignment.tajimas_d(100,3)
snp_map = structure.map(snp_data, ref=ref_sequence, radius=15, force_match=False, output='atom') # If reference seq not defined, then will base it off PDB residue numbering
tajimasd_map = structure.map('tajimas_d', align_file, ref=None, radius=15, force_match=False, output='atom') #Will use first sequence in alignment as a reference sequence.
#Or could use an aligned reference seq.
#Can also accept custom methods. Custom method accepts multiple sequence alignment object and returns a value/object.
#Output can either be by residue or by atom numbering. Residue is default.
#Output is a dictionary mapping residue/atom number to output object.

#OR

mapobject = structmap.Protein(pdbfile,alignfile,refsequence)

tajd = mapobject.tajimas_d(100,3)
tajd_struct = mapobject.tajimas_d_pdb(radius=15)
