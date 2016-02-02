from unittest import TestCase

from structmap import utils, structmap, seqtools
from Bio import AlignIO

class TestUtils(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_is_unambig_dna(self):
        self.assertTrue(structmap.utils.is_unambig_dna("ACTGTCTGTCAT"))
        #Test casing
        self.assertTrue(structmap.utils.is_unambig_dna("actGTCTGtcaT"))
        self.assertFalse(structmap.utils.is_unambig_dna("ACTGMYTGTCAT"))
        self.assertFalse(structmap.utils.is_unambig_dna("ACTG ACTC"))
        self.assertFalse(structmap.utils.is_unambig_dna("ACTG.ACTC"))

    def test_to_string(self):
        pass

class TestPdbtools(TestCase):
    pass

class TestSeqtools(TestCase):
    def setUp(self):
        self.alignment = AlignIO.read('./tests/msa/MSA_test.fsa', 'fasta')
        self.varsites = seqtools._var_site(self.alignment)

    def tearDown(self):
        pass

    def test_var_site(self):
        varsites_keys_to_match = [28, 46, 67, 93, 95, 98, 100]
        self.assertEqual(sorted(self.varsites.keys()),varsites_keys_to_match)

    def test_join_alignments(self):
        msa1 = self.alignment[:,1:2]
        msa2 = self.alignment[:,2:3]
        msa3 = self.alignment[:,3:6]
        d = {}
        for i, msa in enumerate([msa1,msa2,msa3]):
            d[i] = msa
        joined = seqtools._join_alignments(d)
        self.assertEqual(joined.format('fasta'),self.alignment[:,1:6].format('fasta'))
        pass


class TestGentests(TestCase):
    pass

class TestStructmap(TestCase):
    def test_sequence_class(self):
        test_file = "./tests/fasta/test_seq.fsa"
        seq_to_match = ("MKCNISIYFFASFFVLYFAKARNEYDIKENEKFLDVYKEKFNELDKKKYGNVQKTDKKIF"
                       "TFIENKLDILNNSKFNKRWKSYGTPDNIDKNMSLINKHNNEEMFNNNYQSFLSTSSLIKQ"
                       "NKYVPINAVRVSRILSFLDSRIN"
                       )
        seqobj = structmap.Sequence(test_file)
        sequence = seqobj.sequence()
        self.assertEqual(sequence,seq_to_match)
