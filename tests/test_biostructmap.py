from __future__ import absolute_import, division, print_function

import io
from unittest import TestCase
import numpy as np

import Bio.PDB
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from biostructmap import biostructmap, seqtools, gentests, pdbtools
from biostructmap.seqtools import (_sliding_window, _sliding_window_var_sites,
                                   _construct_sub_align)
from biostructmap.gentests import _tajimas_d
from biostructmap.pdbtools import _euclidean_distance_matrix
from biostructmap.map_functions import _tajimas_d

STANDARD_AA_3_LETTERS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                         'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                         'THR', 'TRP', 'TYR', 'VAL']

class TestPdbtools(TestCase):
    def setUp(self):
        self.test_pdb_file = './tests/pdb/1zrl.pdb'
        parser = Bio.PDB.PDBParser()
        pdbname = self.test_pdb_file
        #Get Bio.PDB structure
        self.test_structure = parser.get_structure(pdbname, self.test_pdb_file)
        self.test_model = self.test_structure[0]
        self.test_chain = self.test_model['A']

    def test_euclidean_matrix_calculation(self):
        mat = _euclidean_distance_matrix(self.test_model, selector="all")
        length_to_match = 4908
        self.assertEqual(len(mat[0]), length_to_match)
        self.assertEqual(len(mat[1]), length_to_match)
        #Check that diagonal is all 0
        check_diagonal = (mat[0].diagonal() == np.zeros(length_to_match)).all()
        self.assertTrue(check_diagonal)

    def test_euclidean_matrix_calculation_with_CA(self):
        mat = _euclidean_distance_matrix(self.test_model, selector="CA")
        length_to_match = 583
        self.assertEqual(len(mat[0]), length_to_match)
        self.assertEqual(len(mat[1]), length_to_match)
        #Check that diagonal is all 0
        check_diagonal = (mat[0].diagonal() == np.zeros(length_to_match)).all()
        self.assertTrue(check_diagonal)

    def test_euclidean_matrix_calculation_with_CB(self):
        mat = _euclidean_distance_matrix(self.test_model, selector="CB")
        length_to_match = 583
        self.assertEqual(len(mat[0]), length_to_match)
        self.assertEqual(len(mat[1]), length_to_match)
        #Check that diagonal is all 0
        check_diagonal = (mat[0].diagonal() == np.zeros(length_to_match)).all()
        self.assertTrue(check_diagonal)

    def test_nearby_residues_function_with_CA(self):
        radius = 15
        selector = 'CA'
        test_residue = ('A',(' ', 57, ' '))
        test_residue_to_match = [48,  49,  50,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,
                                 62,  63,  64,  65,  66,  98,  99, 101, 102, 103, 104, 116, 117,
                                 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 136,
                                 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 151]
        residues_to_match = {('A',(' ', x, ' ')) for x in test_residue_to_match}
        nearby = pdbtools.nearby(self.test_model, radius, selector)
        result = nearby[test_residue]
        self.assertEqual(result, residues_to_match)

    def test_nearby_residues_function_with_all_atoms(self):
        radius = 15
        selector = 'all'
        test_residue = ('A',(' ', 57, ' '))
        test_residue_to_match = [8, 11, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
                                 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                                 69, 98, 99, 100, 101, 102, 103, 104, 105, 115, 116, 117, 118,
                                 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
                                 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
                                 145, 146, 147, 148, 149, 150, 151, 152, 154]
        residues_to_match = {('A',(' ', x, ' ')) for x in test_residue_to_match}
        nearby = pdbtools.nearby(self.test_model, radius, selector)
        result = nearby[test_residue]
        self.assertEqual(result, residues_to_match)

    def test_get_pdb_sequence(self):
        filename = './tests/pdb/1zrl.pdb'
        sequence = pdbtools.get_pdb_seq(filename)
        to_match = {'A': (
            'GRQTSSNNEVLSNCREKRKGMKWDCKKKNDRSNYVCIPDRRIQLCIVNLAII'
            'KTYTKETMKDHFIEASKKESQLLLKKNDNKYNSKFCNDLKNSFLDYGHLAMGN'
            'DMDFGGYSTKAENKIQEVFKGAHGEISEHKIKNFRKKWWNEFREKLWEAMLSE'
            'HKNNINNCKNIPQEELQITQWIKEWHGEFLL'
            'ERDNRAKLPKSKCKNNALYEACEKECIDPCMKYRDWIIRSKFEWHTLSKEYETQKVPKENAEN'
            'YLIKISENKNDAKVSLLLNNCDAEYSKYCDCKHTTTLVKSVLNGNDNTIKEKREHIDLDDFSK'
            'FGCDKNSVDTNTKVWECKKPYKLSTKDVCVPPRRQELCLGNIDRIYDKNLLMIKEHILAIAIY'
            'ESRILKRKYKNKDDKEVCKIINKTFADIRDIIGGTDYWNDLSNRKLVGKINTNSNYVHRNKQN'
            'DKLFRDEWWKVIKKDVWNVISWVFKDKTVCKEDDIENIPQFFRWFSEWGDDYCQDKTKMIETL'
            'KVECKEKPCEDDNCKRKCNSYKEWISKKKEEYNKQAKQYQEYQKGNNYKMYSEFKSIKPEVYL'
            'KKYSEKCSNLNFEDEFKEELHSDYKNKCTMCPEVK')
              }
        self.assertEqual(sequence, to_match)
        with self.assertRaises(IOError):
            pdbtools.get_pdb_seq('not_a_file')

    def test_get_mmcif_seq(self):
        filename = './tests/pdb/4nuv.cif'
        mmcif_dict = MMCIF2Dict(filename)
        sequence = pdbtools.get_mmcif_canonical_seq(mmcif_dict)
        to_match = {'C': 'GPTGTENSSQLDFEDVWNSSYGVNDSFPDGDYGA',
                    'D': 'GPTGTENSSQLDFEDVWNSSYGVNDSFPDGDYGA',
                    'A': ('ASNTVMKNCNYKRKRRERDWDCNTKKDVCIPDRRYQLCMKELTNLVNNTDT'
                          'NFHRDITFRKLYLKRKLIYDAAVEGDLLLKLNNYRYNKDFCKDIRWSLGDF'
                          'GDIIMGTDMEGIGYSKVVENNLRSIFGTDEKAQQRRKQWWNESKAQIWTAM'
                          'MYSVKKRLKGNFIWICKLNVAVNIEPQIYRWIREWGRDYVSELPTEVQKLK'
                          'EKCDGKINYTDKKVCKVPPCQNACKSYDQWITRKKNQWDVLSNKFISVKNA'
                          'EKVQTAGIVTPYDILKQELDEFNEVAFENEINKRDGAYIELCVCSVEEAKK'
                          'NTQEVVTNVDN'),
                    'B':  ('ASNTVMKNCNYKRKRRERDWDCNTKKDVCIPDRRYQLCMKELTNLVNNTDT'
                           'NFHRDITFRKLYLKRKLIYDAAVEGDLLLKLNNYRYNKDFCKDIRWSLGDF'
                           'GDIIMGTDMEGIGYSKVVENNLRSIFGTDEKAQQRRKQWWNESKAQIWTAM'
                           'MYSVKKRLKGNFIWICKLNVAVNIEPQIYRWIREWGRDYVSELPTEVQKLK'
                           'EKCDGKINYTDKKVCKVPPCQNACKSYDQWITRKKNQWDVLSNKFISVKNA'
                           'EKVQTAGIVTPYDILKQELDEFNEVAFENEINKRDGAYIELCVCSVEEAKK'
                           'NTQEVVTNVDN')}
        self.assertDictEqual(sequence, to_match)

    def test_tajimas_d_on_structure(self):
        #test_sequence_alignment = AlignIO.read('./tests/msa/msa_test_86-104', 'fasta')
        test_sequence_alignment = {('A',): biostructmap.SequenceAlignment(
            './tests/msa/msa_test_86-104', 'fasta')}
        test_ref_dict = {('A', (' ', x+86, ' ')): ('A', (x*3+1, x*3 + 2, x*3 + 3)) for
                         x in range(0, 18)}
        test_surrounding_residues = [('A', (' ', x, ' ')) for x in range(86, 104)]
        result = _tajimas_d(self.test_structure, test_sequence_alignment,
                            test_surrounding_residues, test_ref_dict)
        self.assertEqual(result, -0.7801229937910628)

    def test_tajimas_d_on_structure_with_subset_of_reference_residues(self):
        #test_sequence_alignment = AlignIO.read('./tests/msa/msa_test_86-104', 'fasta')
        test_sequence_alignment = {('A',): biostructmap.SequenceAlignment(
            './tests/msa/msa_test_86-104', 'fasta')}
        test_ref_dict = {('A', (' ', x+86, ' ')): ('A', (x*3 + 1, x*3 + 2, x*3 + 3)) for
                         x in range(18)}
        test_surrounding_residues = [('A', (' ', x, ' ')) for x in range(86, 96)]
        result = _tajimas_d(self.test_structure, test_sequence_alignment,
                            test_surrounding_residues, test_ref_dict)
        self.assertEqual(result, -0.709896167879475)

class TestSeqtools(TestCase):
    def setUp(self):
        self.test_file = './tests/msa/MSA_test.fsa'
        self.alignment = AlignIO.read(self.test_file, 'fasta')
        self.biostructmap_alignment = biostructmap.SequenceAlignment('./tests/msa/MSA_test.fsa')
        self.varsites = seqtools._var_site(self.alignment)

    def tearDown(self):
        pass

    def test_var_site(self):
        varsites_keys_to_match = [28, 46, 67, 93, 95, 98, 100]
        self.assertEqual(sorted(self.varsites.keys()), varsites_keys_to_match)

    def test_join_alignments(self):
        msa1 = self.alignment[:, 1:2]
        msa2 = self.alignment[:, 2:3]
        msa3 = self.alignment[:, 3:6]
        alignment_dict = {}
        for i, msa in enumerate([msa1, msa2, msa3]):
            alignment_dict[i] = msa
        joined = seqtools._join_alignments(alignment_dict)
        self.assertEqual(joined.format('fasta'),
                         self.alignment[:, 1:6].format('fasta'))

    def test_sliding_window(self):
        step = 3
        length = 10
        #Test sliding window output in AlignIO format
        slider = _sliding_window(self.alignment, length, step=step, fasta_out=False)
        for i, window in enumerate(slider):
            to_equal = self.alignment[:, i*step:i*step+length]
            to_equal = to_equal.format('fasta')
            window_fasta = window.format('fasta')
            self.assertEqual(window_fasta, to_equal)
        #Test fasta output
        slider = _sliding_window(self.alignment, length, step=step, fasta_out=True)
        for i, window in enumerate(slider):
            to_equal = self.alignment[:, i*step:i*step+length]
            to_equal = to_equal.format('fasta')
            self.assertEqual(window, to_equal)
        #Test reading from file
        slider = _sliding_window(self.test_file, length, step=step, fasta_out=False)
        for i, window in enumerate(slider):
            to_equal = self.alignment[:, i*step:i*step+length]
            to_equal = to_equal.format('fasta')
            window_fasta = window.format('fasta')
            self.assertEqual(window_fasta, to_equal)
        #Test with a different step size
        slider = _sliding_window(self.alignment, length, step=1, fasta_out=False)
        for i, window in enumerate(slider):
            to_equal = self.alignment[:, i:i+length]
            to_equal = to_equal.format('fasta')
            window_fasta = window.format('fasta')
            self.assertEqual(window_fasta, to_equal)

    def test_sliding_window_var_sites(self):
        step = 3
        length = 10
        slider = _sliding_window_var_sites(self.alignment, length, step)
        varsites_keys_to_match = [28, 46, 67, 93, 95, 98, 100]
        null_align = self.alignment[:, 0:0]
        for i, window in enumerate(slider):
            #Check if key within range
            in_range = [x for x in varsites_keys_to_match if (step*i) <= x < (step*i + length)]
            print(in_range)
            if in_range:
                window_i = self.alignment[:, in_range[0]:in_range[0]+1]
                if len(in_range) > 1:
                    for x in in_range[1:]:
                        window_i = window_i + self.alignment[:, x:x+1]
                self.assertEqual(window_i.format('fasta'),
                                 window.format('fasta'))
            else:
                self.assertEqual(window.format('fasta'),
                                 null_align.format('fasta'))

    def test_sliding_window_var_sites_with_file(self):
        step = 3
        length = 10
        slider = _sliding_window_var_sites(self.test_file, length, step)
        varsites_keys_to_match = [28, 46, 67, 93, 95, 98, 100]
        null_align = self.alignment[:, 0:0]
        for i, window in enumerate(slider):
            #Check if key within range
            in_range = [x for x in varsites_keys_to_match if (step*i) <= x < (step*i + length)]
            print(in_range)
            if in_range:
                window_i = self.alignment[:, in_range[0]:in_range[0]+1]
                if len(in_range) > 1:
                    for x in in_range[1:]:
                        window_i = window_i + self.alignment[:, x:x+1]
                self.assertEqual(window_i.format('fasta'),
                                 window.format('fasta'))
            else:
                self.assertEqual(window.format('fasta'),
                                 null_align.format('fasta'))

    def test_blast_sequences(self):
        seq1 = "GSNAKFGLWVDGNCEDIPHVNEFPAID"
        seq1_bio = Seq(seq1)
        seq2 = "NAKFGLWV"
        seq2_bio = Seq(seq2)
        test_map_forward, test_map_reverse = seqtools.blast_sequences(seq1, seq2)
        forward_match = {3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8}
        reverse_match = {1: 3, 2: 4, 3: 5, 4: 6, 5: 7, 6: 8, 7: 9, 8: 10}
        #Check alignment for string input
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)
        #Check alignment for Bio.Seq input
        test_map_forward, test_map_reverse = seqtools.blast_sequences(seq1_bio, seq2_bio)
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)

    def test_pairwise_align_sequences(self):
        seq1 = "GSNAKFGLWVDGNCEDIPHVNEFPAID"
        seq1_bio = Seq(seq1)
        seq2 = "NAKFGLWV"
        seq2_bio = Seq(seq2)
        test_map_forward, test_map_reverse = seqtools.pairwise_align(seq1, seq2)
        forward_match = {3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8}
        reverse_match = {1: 3, 2: 4, 3: 5, 4: 6, 5: 7, 6: 8, 7: 9, 8: 10}
        #Check alignment for string input
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)
        #Check alignment for Bio.Seq input
        test_map_forward, test_map_reverse = seqtools.pairwise_align(seq1_bio, seq2_bio)
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)

    def test_blast_sequences_with_gaps_in_first_sequence(self):
        '''Test with reverse input (ie gaps in first sequence)'''
        seq1 = "GSNAKFGLWVDGNCEDIPHVNEFPAID"
        seq1_bio = Seq(seq1)
        seq2 = "NAKFLWVDG"
        seq2_bio = Seq(seq2)
        test_map_reverse, test_map_forward = seqtools.blast_sequences(seq2, seq1)
        forward_match = {3: 1, 4: 2, 5: 3, 6: 4, 8: 5, 9: 6, 10: 7, 11: 8, 12: 9}
        reverse_match = {1: 3, 2: 4, 3: 5, 4: 6, 5: 8, 6: 9, 7: 10, 8: 11, 9: 12}
        #Check alignment for string input
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)
        #Check alignment for Bio.Seq input
        test_map_forward, test_map_reverse = seqtools.blast_sequences(seq1_bio, seq2_bio)
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)

    def test_pairwise_align_sequences_with_gaps_in_first_sequence(self):
        '''Test with reverse input (ie gaps in first sequence)'''
        seq1 = "GSNAKFGLWVDGNCEDIPHVNEFPAID"
        seq1_bio = Seq(seq1)
        seq2 = "NAKFLWVDG"
        seq2_bio = Seq(seq2)
        test_map_reverse, test_map_forward = seqtools.pairwise_align(seq2, seq1)
        forward_match = {3: 1, 4: 2, 5: 3, 6: 4, 8: 5, 9: 6, 10: 7, 11: 8, 12: 9}
        reverse_match = {1: 3, 2: 4, 3: 5, 4: 6, 5: 8, 6: 9, 7: 10, 8: 11, 9: 12}
        #Check alignment for string input
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)
        #Check alignment for Bio.Seq input
        test_map_forward, test_map_reverse = seqtools.pairwise_align(seq1_bio, seq2_bio)
        self.assertEqual(forward_match, test_map_forward)
        self.assertEqual(reverse_match, test_map_reverse)

    def test_protein_dna_alignment(self):
        '''Need to account for an intron frameshift when converting from DNA to
        protein sequence.
        '''
        dna = 'ATGAAATGTAATATTAGTATATATTTTTTTATGAAATGTAATATTAGTATATATTTTTTT'
        protein = 'MKCNISIYFFMKCNISIYFF'
        result = seqtools.align_protein_to_dna(protein, dna)
        result_keys_to_match = [i for i in range(1, 21)]
        result_codons_to_match = [(i*3-2, i*3-1, i*3) for i in range(1, 21)]
        self.assertEqual(result_keys_to_match, sorted(result))
        sorted_codons = [result[i] for i in sorted(result)]
        self.assertEqual(result_codons_to_match, sorted_codons)

    def test_protein_dna_alignment_without_exonerate(self):
        '''Need to account for an intron frameshift when converting from DNA to
        protein sequence.
        '''
        seqtools.LOCAL_EXONERATE = False
        dna = 'ATGAAATGTAATATTAGTATATATTTTTTTATGAAATGTAATATTAGTATATATTTTTTT'
        protein = 'MKCNISIYFFMKCNISIYFF'
        result = seqtools.align_protein_to_dna(protein, dna)
        result_keys_to_match = [i for i in range(1, 21)]
        result_codons_to_match = [(i*3-2, i*3-1, i*3) for i in range(1, 21)]
        self.assertEqual(result_keys_to_match, sorted(result))
        sorted_codons = [result[i] for i in sorted(result)]
        self.assertEqual(result_codons_to_match, sorted_codons)
        seqtools.LOCAL_EXONERATE = True

    def test_protein_dna_alignment_with_intron(self):
        '''Need to account for an intron frameshift when converting from DNA to
        protein sequence.
        '''
        dna = ('ATGAAATGTAATATTAGTATATATTTTTTT'
               'GTTGTATAG'
               'ATGAAATGTAATATTAGTATATATTTTTTTATGAAATGTAATATTAGTATATATTTTTTT')
        protein = 'MKCNISIYFFMKCNISIYFFMKCNISIYFF'
        result = seqtools.align_protein_to_dna(protein, dna)
        result_keys_to_match = [i for i in range(1,11)] + [i for i in range(11, 31)]
        result_codons_to_match = [(i*3-2, i*3-1, i*3) for i in range(1, 11)] + \
                [(i*3-2, i*3-1, i*3) for i in range(14, 34)]
        self.assertEqual(result_keys_to_match, sorted(result))
        sorted_codons = [result[i] for i in sorted(result)]
        self.assertEqual(result_codons_to_match, sorted_codons)

    def test_sub_align(self):
        codons = [(1, 2, 3), (4, 5, 6), (10, 11, 12)]
        result = _construct_sub_align(self.biostructmap_alignment, codons, 'fasta')
        to_match = self.alignment[:, 0:6] + self.alignment[:, 9:12]
        self.assertEqual(to_match.format('fasta'), result)

class TestGentests(TestCase):
    def setUp(self):
        self.test_file = './tests/msa/MSA_test.fsa'
        self.small_test_file = './tests/msa/MSA_small.fsa'
        self.alignment = AlignIO.read(self.test_file, 'fasta')
        self.small_alignment = AlignIO.read(self.small_test_file, 'fasta')
    def test_tajimas_d(self):
        #Test basic calculation of Tajima's D
        taj_d = gentests.tajimas_d(self.alignment)
        self.assertEqual(taj_d, -1.553110875316991)
        #Test Tajima's D over a sliding window
        to_match = [(25.5, -1.2371599802089934),
                    (37.5, -1.2371599802089934),
                    (49.5, -1.3584148210238671),
                    (61.5, -1.2371599802089934),
                    (73.5, -1.3584148210238671),
                    (85.5, -1.486139836697004),
                    (97.5, -1.434138485438007),
                    (109.5, -1.434138485438007),
                    (121.5, -1.2371599802089934)]
        taj_d_window = gentests.tajimas_d(self.alignment, 50, 12)
        taj_d_window = [(x, taj_d_window[x]) for x in sorted(taj_d_window)]
        self.assertEqual(taj_d_window, to_match)
        #Test that method returns 'None' when Tajimas D can't be calculated.
        taj_d_small_window = gentests.tajimas_d(self.alignment, 12, 3)
        self.assertEqual(taj_d_small_window[6.5], None)

    def test_tajimas_d_with_divide_by_zero_error(self):
        taj_d = gentests._tajimas_d(self.small_alignment)
        self.assertEqual(taj_d, None)

    def test_nucleotide_diversity_calculation(self):
        diversity = gentests.nucleotide_diversity(self.small_alignment)
        self.assertAlmostEqual(diversity, 2 / (3 * 31))

    def test_wattersons_theta_calculation(self):
        theta = gentests.wattersons_theta(self.small_alignment)
        self.assertAlmostEqual(theta, 1/1.5)

    def test_shannon_entropy(self):
        entropy = gentests.shannon_entropy(self.small_alignment[:,0:27])
        self.assertAlmostEqual(entropy, 0.10203287045049884)

    def test_normalized_shannon_entropy(self):
        entropy = gentests.shannon_entropy(self.small_alignment[:,0:27], normalized=True)
        self.assertAlmostEqual(entropy, 0.023608183248397613)

    def test_shannon_entropy_basic_calculation(self):
        entropy = gentests._calculate_shannon_entropy(['A', 'B', 'C', 'D'],
                                                      protein_letters='ABCD')
        self.assertAlmostEqual(entropy, 2)
        entropy = gentests._calculate_shannon_entropy(['A', 'A', 'A', 'A'],
                                                      protein_letters='ABCD')
        self.assertAlmostEqual(entropy, 0)

class Testbiostructmap(TestCase):
    def setUp(self):
        self.test_file = './tests/pdb/1as5.pdb'
        self.test_align = './tests/msa/MSA_small.fsa'

    def test_structure_object_instantiation(self):
        structure = biostructmap.Structure(self.test_file)
        for model in structure:
            #Test iteration
            self.assertTrue(isinstance(model, biostructmap.Model))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(structure[model.get_id()], biostructmap.Model))
        self.assertTrue(isinstance(structure.structure, Bio.PDB.Structure.Structure))
        self.assertTrue(isinstance(structure.structure[0]['A'], Bio.PDB.Chain.Chain))
        self.assertTrue(isinstance(structure.sequences, dict))

    def test_structure_object_instantiation_with_file_like_object(self):
        with open(self.test_file, 'r') as f:
            test_filelike = io.StringIO(f.read())
        structure = biostructmap.Structure(test_filelike)
        for model in structure:
            #Test iteration
            self.assertTrue(isinstance(model, biostructmap.Model))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(structure[model.get_id()], biostructmap.Model))
        self.assertTrue(isinstance(structure.structure, Bio.PDB.Structure.Structure))
        self.assertTrue(isinstance(structure.sequences, dict))


    def test_model_instantiation(self):
        structure = biostructmap.Structure(self.test_file)
        test_model = structure[0]
        self.assertEqual(test_model.get_id(),0)
        self.assertTrue(isinstance(test_model.model, Bio.PDB.Model.Model))
        self.assertEqual(test_model.parent(), structure)
        for chain in test_model:
            #Test iteration
            self.assertTrue(isinstance(chain, biostructmap.Chain))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(test_model[chain.get_id()], biostructmap.Chain))

    def test_chain_instantiation(self):
        structure = biostructmap.Structure(self.test_file)
        test_chain = structure[0]['A']
        self.assertEqual(test_chain.get_id(),'A')
        self.assertTrue(isinstance(test_chain.chain, Bio.PDB.Chain.Chain))
        self.assertEqual(test_chain.parent(), structure[0])
        for residue in test_chain:
            #Test iteration
            self.assertTrue(isinstance(residue, Bio.PDB.Residue.Residue))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(test_chain[residue.get_id()], Bio.PDB.Residue.Residue))

    def test_structure_nearby_function(self):
        structure = biostructmap.Structure(self.test_file)
        result = structure.nearby()
        self.assertTrue(isinstance(result, dict))
        for i in result.values():
            self.assertTrue(isinstance(i, set))

    def test_rsa_determination(self):
        chain = biostructmap.Structure(self.test_file)[0]['A']
        result = chain.rel_solvent_access()
        for residue in [residues for residues in chain if residues.get_id()[0] == ' ']:
            self.assertTrue(isinstance(result[residue.get_id()], float))

    def test_default_mapping_procedure(self):
        structure = biostructmap.Structure(self.test_file)
        chain = biostructmap.Structure(self.test_file)[0]['A']
        data = {'A': [x for x in range(0, 25)]}
        mapping = structure.map(data)
        for residue in [residues for residues in chain if residues.get_id()[0] == ' ']:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertTrue(isinstance(result, float))
        mapping = structure.map(data, method='default', ref=None, radius=0, selector='all')
        nondetermined_residues = [1, 2, 3]
        for residue in [residues for residues in chain if
                        residues.get_id()[0] == ' ' and
                        residues.get_id()[1] not in nondetermined_residues]:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertEqual(result, residue.get_id()[1])
        #Test if data is in a dictionary
        mapping = structure.map({'A': {x:x for x in range(0, 25)}})
        for residue in [residues for residues in chain if
                        residues.get_id()[0] == ' ' and
                        residues.get_id()[1] not in nondetermined_residues]:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertTrue(isinstance(result, float))
        #Test that an rsa_range of 0-1 doesn't change results
        rsa_mapping = structure.map(data, rsa_range=[0,1])
        self.assertDictEqual(rsa_mapping, mapping)

    def test_default_mapping_procedure_with_pairwise(self):
        seqtools.LOCAL_BLAST = False
        structure = biostructmap.Structure(self.test_file)
        chain = biostructmap.Structure(self.test_file)[0]['A']
        data = {'A': [x for x in range(0, 25)]}
        mapping = structure.map(data)
        for residue in [residues for residues in chain if residues.get_id()[0] == ' ']:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertTrue(isinstance(result, float))
        mapping = structure.map(data, method='default', ref=None, radius=0, selector='all')
        nondetermined_residues = [1, 2, 3]
        for residue in [residues for residues in chain if
                        residues.get_id()[0] == ' ' and
                        residues.get_id()[1] not in nondetermined_residues]:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertEqual(result, residue.get_id()[1])
        #Test if data is in a dictionary
        mapping = structure.map({'A': {x:x for x in range(0, 25)}})
        for residue in [residues for residues in chain if
                        residues.get_id()[0] == ' ' and
                        residues.get_id()[1] not in nondetermined_residues]:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertTrue(isinstance(result, float))
        #Test that an rsa_range of 0-1 doesn't change results
        rsa_mapping = structure.map(data, rsa_range=[0,1])
        self.assertDictEqual(rsa_mapping, mapping)
        seqtools.LOCAL_BLAST = True

    def test_rsa_filtering_procedure(self):
        data = {'A': [x for x in range(0, 25)]}
        structure = biostructmap.Structure(self.test_file)
        mapping = structure.map(data)
        filtered_mapping = structure.map(data, rsa_range=[0.2, 1])
        self.assertEqual(filtered_mapping[('A', (' ', 13, ' '))], None)
        self.assertEqual(set(mapping.keys()) - set([x for x in
                                                    filtered_mapping.keys() if
                                                    filtered_mapping[x] is not None]),
                         set([('A', (' ', 4 , ' ')), ('A', (' ', 13, ' ')),
                              ('A', (' ', 16, ' '))]))


    def test_sequence_alignment_instantiation(self):
        test_align = biostructmap.SequenceAlignment(self.test_align)
        self.assertTrue(isinstance(test_align.alignment,
                                   Bio.Align.MultipleSeqAlignment))
        for i, seq in enumerate(test_align):
            #Test iteration
            self.assertTrue(isinstance(seq, Bio.SeqRecord.SeqRecord))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(test_align[i], Bio.SeqRecord.SeqRecord))

    def test_tajimas_d_on_sequence_alignment(self):
        test_align = biostructmap.SequenceAlignment('./tests/msa/MSA_test.fsa')
        #Test basic calculation of Tajima's D
        taj_d = test_align.tajimas_d()
        self.assertEqual(taj_d, -1.553110875316991)
        #Test for errors in input to the function
        with self.assertRaises(TypeError):
            test_align.tajimas_d(3.5,5)
        with self.assertRaises(TypeError):
            test_align.tajimas_d(3,5.5)

    def test_tajimas_d_on_long_sequence(self):
        test_align = biostructmap.SequenceAlignment('./tests/msa/MSA_test_long.fsa')
        taj_d = test_align.tajimas_d()
        self.assertEqual(taj_d, 0.33458440732186856)

    def test_residue_number_to_atom_number_function(self):
        structure = biostructmap.Structure(self.test_file)
        mapping = structure[0]['A'].residue_to_atom_map()
        self.assertEqual(len(mapping), 25)
        first_residue_atoms = list(range(1,21))
        self.assertEqual(sorted(mapping[1]), first_residue_atoms)

    def test_secondary_structure_dictionary_creation(self):
        chain = biostructmap.Structure(self.test_file)[0]['A']
        ss_dict_numeric = {1: 7, 2: 7, 3: 7, 4: 6, 5: 6, 6: 6,
                           7: 5, 8: 5, 9: 7, 10: 7, 11: 7, 12: 7,
                           13: 7, 14: 5, 15: 5, 16: 7, 17: 5, 18: 5,
                           19: 7, 20: 7, 21: 7, 22: 6, 23: 7, 24: 7}
        ss_dict = {1: '-', 2: '-', 3: '-', 4: 'S', 5: 'S',
                   6: 'S', 7: 'T', 8: 'T', 9: '-', 10: '-',
                   11: '-', 12: '-', 13: '-', 14: 'T', 15: 'T',
                   16: '-', 17: 'T', 18: 'T', 19: '-', 20: '-',
                   21: '-', 22: 'S', 23: '-', 24: '-'}
        _ss_dict_numeric = {(' ', x, ' '): y for x, y in ss_dict_numeric.items()}
        _ss_dict = {(' ', x, ' '): y for x, y in ss_dict.items()}
        self.assertDictEqual(_ss_dict, chain.secondary_structure())
        self.assertDictEqual(_ss_dict_numeric,
                             chain.secondary_structure(numeric_ss_code=True))


class TestbiostructmapMMCIF(TestCase):
    def setUp(self):
        self.test_file = './tests/pdb/1as5.cif'
        self.test_align = './tests/msa/MSA_small.fsa'
        self.structure = biostructmap.Structure(self.test_file, mmcif=True)

    def test_structure_object_instantiation(self):
        structure = self.structure
        for model in structure:
            #Test iteration
            self.assertTrue(isinstance(model, biostructmap.Model))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(structure[model.get_id()], biostructmap.Model))
        self.assertTrue(isinstance(structure.structure, Bio.PDB.Structure.Structure))
        self.assertTrue(isinstance(structure.structure[0]['A'], Bio.PDB.Chain.Chain))
        self.assertTrue(isinstance(structure.sequences, dict))

    def test_structure_object_instantiation_with_file_like_object(self):
        with open(self.test_file, 'r') as f:
            test_filelike = io.StringIO(f.read())
        structure = biostructmap.Structure(test_filelike, mmcif=True)
        for model in structure:
            #Test iteration
            self.assertTrue(isinstance(model, biostructmap.Model))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(structure[model.get_id()], biostructmap.Model))
        self.assertTrue(isinstance(structure.structure, Bio.PDB.Structure.Structure))
        self.assertTrue(isinstance(structure.sequences, dict))


    def test_model_instantiation(self):
        structure = self.structure
        test_model = structure[0]
        self.assertEqual(test_model.get_id(),0)
        self.assertTrue(isinstance(test_model.model, Bio.PDB.Model.Model))
        self.assertEqual(test_model.parent(), structure)
        for chain in test_model:
            #Test iteration
            self.assertTrue(isinstance(chain, biostructmap.Chain))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(test_model[chain.get_id()], biostructmap.Chain))

    def test_chain_instantiation(self):
        structure = self.structure
        test_chain = structure[0]['A']
        self.assertEqual(test_chain.get_id(), 'A')
        self.assertTrue(isinstance(test_chain.chain, Bio.PDB.Chain.Chain))
        self.assertEqual(test_chain.parent(), structure[0])
        for residue in test_chain:
            #Test iteration
            self.assertTrue(isinstance(residue, Bio.PDB.Residue.Residue))
            #test _getitem__ method work to return a model object
            self.assertTrue(isinstance(test_chain[residue.get_id()], Bio.PDB.Residue.Residue))

    def test_structure_nearby_function(self):
        structure = self.structure
        result = structure.nearby()
        self.assertTrue(isinstance(result, dict))
        for i in result.values():
            self.assertTrue(isinstance(i, set))

    def test_rsa_determination(self):
        chain = self.structure[0]['A']
        result = chain.rel_solvent_access()
        for residue in chain:
            if residue.get_resname() not in STANDARD_AA_3_LETTERS:
                self.assertEqual(result.get(residue.get_id(), None), None)
            else:
                self.assertTrue(isinstance(result[residue.get_id()], float))

    def test_default_mapping_procedure(self):
        structure = self.structure
        chain = self.structure[0]['A']
        data = {'A': [x for x in range(0, 25)]}
        mapping = structure.map(data)
        for residue in [residues for residues in chain if residues.get_id()[0] == ' ']:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertTrue(isinstance(result, float))
        mapping = structure.map(data, method='default', ref=None, radius=0, selector='all')
        nondetermined_residues = [25]
        for residue in [residues for residues in chain if
                        residues.get_id()[0] == ' ' and
                        residues.get_id()[1] not in nondetermined_residues]:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertEqual(result, residue.get_id()[1])
        #Test if data is in a dictionary
        mapping = structure.map({'A': {x:x for x in range(0, 25)}})
        for residue in [residues for residues in chain if
                        residues.get_id()[0] == ' ' and
                        residues.get_id()[1] not in nondetermined_residues]:
            result = mapping[residue.get_full_id()[2:4]]
            self.assertTrue(isinstance(result, float))

    def test_rsa_filtering_procedure(self):
        data = {'A': [x for x in range(0, 25)]}
        structure = self.structure
        mapping = structure.map(data)
        filtered_mapping = structure.map(data, rsa_range=[0.2, 1])
        self.assertEqual(filtered_mapping[('A', (' ', 13, ' '))], None)
        self.assertEqual(set(mapping.keys()) - set([x for x in
                                                    filtered_mapping.keys() if
                                                    filtered_mapping[x] is not None]),
                         set([('A', (' ', 4, ' ')), ('A', (' ', 13, ' ')),
                              ('A', (' ', 25, ' ')), ('A', (' ', 14, ' ')),
                              ('A', (' ', 3, ' ')), ('A', (' ', 2, ' ')),
                              ('A', (' ', 16, ' '))]))


    def test_residue_number_to_atom_number_function(self):
        '''
        Residue to atom map won't function with current Biopython version
        '''
        structure = self.structure
        with self.assertRaises(TypeError):
            mapping = structure[0]['A'].residue_to_atom_map()

    def test_secondary_structure_dictionary_creation(self):
        chain = self.structure[0]['A']
        ss_dict_numeric = {1: 7, 2: 7, 3: 7, 4: 6, 5: 6, 6: 6,
                           7: 5, 8: 5, 9: 7, 10: 7, 11: 7, 12: 7,
                           13: 7, 14: 5, 15: 5, 16: 7, 17: 5, 18: 5,
                           19: 7, 20: 7, 21: 7, 22: 6, 23: 7, 24: 7}
        ss_dict = {1: '-', 2: '-', 3: '-', 4: 'S', 5: 'S',
                   6: 'S', 7: 'T', 8: 'T', 9: '-', 10: '-',
                   11: '-', 12: '-', 13: '-', 14: 'T', 15: 'T',
                   16: '-', 17: 'T', 18: 'T', 19: '-', 20: '-',
                   21: '-', 22: 'S', 23: '-', 24: '-'}
        _ss_dict_numeric = {(' ', x, ' '): y for x, y in ss_dict_numeric.items()}
        _ss_dict = {(' ', x, ' '): y for x, y in ss_dict.items()}
        self.assertDictEqual(_ss_dict, chain.secondary_structure())
        self.assertDictEqual(_ss_dict_numeric,
                             chain.secondary_structure(numeric_ss_code=True))


class TestStructureMapMethod(TestCase):
    def setUp(self):
        self.test_file = './tests/pdb/4nuv.cif'
        self.structure = biostructmap.Structure(self.test_file, mmcif=True)

    def test_default_mapping_on_multchain_structure(self):
        ref_seq = ('ASNTVMKNCNYKRKRRERDWDCNTKKDVCIPDRRYQLCMKELTNLVNNTDT'
                   'NFHRDITFRKLYLKRKLIYDAAVEGDLLLKLNNYRYNKDFCKDIRWSLGDF'
                   'GDIIMGTDMEGIGYSKVVENNLRSIFGTDEKAQQRRKQWWNESKAQIWTAM'
                   'MYSVKKRLKGNFIWICKLNVAVNIEPQIYRWIREWGRDYVSELPTEVQKLK'
                   'EKCDGKINYTDKKVCKVPPCQNACKSYDQWITRKKNQWDVLSNKFISVKNA'
                   'EKVQTAGIVTPYDILKQELDEFNEVAFENEINKRDGAYIELCVCSVEEAKK'
                   'NTQEVVTNVDN')
        reference_seqs = {'A': ref_seq, 'B': ref_seq}
        data = {('A', 'B'): range(500)}

        mapped = self.structure.map(data=None, method='count_residues',
                                    ref=reference_seqs, radius=3)
        # Number of surrounding residues verified using selections in Pymol
        self.assertEqual(mapped[('A', (' ', 271, ' '))], 10)

        mapped = self.structure.map(data=data, method='snps',
                                    ref=reference_seqs, radius=3,
                                    method_params={'ignore_duplicates':True,
                                                   'output_count': True})
        self.assertEqual(mapped[('A', (' ', 271, ' '))], 9)
        mapped = self.structure.map(data=data, method='snps',
                                    ref=reference_seqs, radius=3,
                                    method_params={'ignore_duplicates':False,
                                                   'output_count': True})
        self.assertEqual(mapped[('A', (' ', 271, ' '))], 10)
