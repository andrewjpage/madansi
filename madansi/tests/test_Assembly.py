import unittest
from madansi.Assembly import Assembly

class TestAssembly(unittest.TestCase):
    def test_empty_file(self):
        assembly = Assembly('madansi/tests/data/empty_file.fa')
        self.assertEqual(assembly.sequence_names(),[])
    
    def test_fasta_file_one_sequence(self):
        assembly = Assembly('madansi/tests/data/fasta_file_one_sequence.fa')
        self.assertEqual(assembly.sequence_names(), ['Sequence1'])
    
    def test_fasta_file_three_sequences(self):
        assembly = Assembly('madansi/tests/data/fasta_file_three_sequences.fa')
        self.assertCountEqual(assembly.sequence_names(), ['Sequence1', 'Sequence3', 'Sequence2'])