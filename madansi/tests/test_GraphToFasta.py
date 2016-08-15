import unittest
from madansi.GraphToFasta import GraphToFasta
import os
import filecmp

class TestGraphToFasta(unittest.TestCase):
    
    def test_empty_graph(self):
        """Should return an empty fasta file"""
        graph = nx.Graph()
        fname = "output.fa"
        graph_to_fasta_object = GraphToFasta(graph, fname)
        self.assertTrue(filecmp.cmp(graph_to_fasta_object.output_file(), 'madansi/tests/data/empty_file.fa'))
        os.unlink(fname)