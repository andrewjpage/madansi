#Tests to check that walk graphs is working properly

import unittest
from madansi.WalkGraphs import WalkGraphs
import filecmp
import networkx as nx
import os


class TestWalkGraphsCases(unittest.TestCase):
    
    def test_one_sequence(self):
        """Tests one sequence with no genes- should return the sequence"""
        output_file = 'output.dot'
        expected_graph = 'madansi/tests/data/WalkGraphs/expected_one_sequence.dot'
        input_file = 'madansi/tests/data/WalkGraphs/one_sequence.dot'
        test_data = 'madansi/tests/data/WalkGraphs/WalkGraphs_test_data'
        wg = WalkGraphs(input_file ,test_data , output_file)
        output_list = wg.create_linear_subgraph()
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(output_file))
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph))
        self.assertTrue(output_list == [])
        self.assertTrue(nx.is_isomorphic(h,g))
        os.unlink(output_file)
        
    def test_two_sequences_no_genes(self):
        """Should return an empty graph with both sequences in the unused sequences list"""
        output_file = 'output.dot'
        expected_graph = 'madansi/tests/data/WalkGraphs/expected_two_sequences_no_genes.dot'
        input_file= 'madansi/tests/data/WalkGraphs/two_sequences_no_genes.dot'
        test_data = 'madansi/tests/data/WalkGraphs/WalkGraphs_test_data'
        wg = WalkGraphs(input_file , test_data, output_file)
        output_list = wg.create_linear_subgraph()
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(output_file))
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph))
        self.assertTrue(output_list == ['Sequence1', 'Sequence2'])
        self.assertTrue(nx.is_isomorphic(h,g))
        os.unlink(output_file)
    
    def test_two_sequences_one_gene(self):
        """Should return a graph with just the sequence with the whole of the sequence with a gene present and the other sequence should go into the unused sequences list"""
        output_file = 'output.dot'
        expected_graph = 'madansi/tests/data/WalkGraphs/expected_two_sequences_one_gene.dot'
        input_file= 'madansi/tests/data/WalkGraphs/two_sequences_one_gene.dot'
        test_data = 'madansi/tests/data/WalkGraphs/WalkGraphs_test_data'
        wg = WalkGraphs(input_file , test_data, output_file)
        output_list = wg.create_linear_subgraph()
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(output_file))
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph))
        self.assertTrue(output_list == ['Sequence2'])
        self.assertTrue(nx.is_isomorphic(h,g))
        os.unlink(output_file)