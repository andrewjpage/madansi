import unittest
from madansi.ContigGraph import ContigGraph
import os
import networkx as nx
import filecmp

class TestContigGraph(unittest.TestCase):
	
    def test__init__(self):
        my_list = []
        my_contig_graph = ContigGraph(my_list)
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/empty_graph.dot'))
        self.assertTrue(nx.is_isomorphic(my_contig_graph.contig_graph, expected_graph))
        
    def test_empty_list(self):
        my_list = []
        my_contig_graph = ContigGraph(my_list)
        my_contig_graph.create_contig_subgraph()
        my_graph = my_contig_graph.contig_graph
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/empty_graph.dot'))
        self.assertTrue(nx.is_isomorphic(my_graph, expected_graph))
    
    def test_one_entry(self):
        my_list = [[('Contig1', 'Contig2'), 3, []]]
        my_contig_graph = ContigGraph(my_list)
        my_contig_graph.create_contig_subgraph()
        my_graph = my_contig_graph.contig_graph
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/contig_graph_one_entry.dot'))
        self.assertTrue(nx.is_isomorphic(my_graph, expected_graph))
        
    def test_three_entries(self):
        my_list = [[('Contig1', 'Contig2'),1, []], [('Contig2', 'Contig3'),2,[]],[ ('Contig3', 'Contig4'),2, []]]
        my_contig_graph = ContigGraph(my_list)
        my_contig_graph.create_contig_subgraph()
        my_graph = my_contig_graph.contig_graph
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/contig_graph_three_entries.dot'))
        self.assertTrue(nx.is_isomorphic(my_graph, expected_graph))
       
    
    def test_five_entries(self):
        my_list = [[('Contig1', 'Contig2'), 2, []], [('Contig3', 'Contig1'),3,[]], [('Contig3','Contig4'),1,[]], [('Contig4', 'Contig1'), 4,[]], [('Contig2', 'Contig5'), 1,[]]]
        my_contig_graph = ContigGraph(my_list)
        my_contig_graph.create_contig_subgraph()
        my_graph = my_contig_graph.contig_graph
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/contig_graph_five_entries.dot'))
        self.assertTrue(nx.is_isomorphic(my_graph, expected_graph))
    
        