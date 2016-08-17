import unittest
import networkx as nx
from madansi.RemoveContigs import RemoveContigs

class TestRemoveContigs(unittest.TestCase):
    
    def test_remove_no_contigs(self):
        my_graph = nx.Graph()
        my_graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        remove_contigs = RemoveContigs(my_graph)
        remove_contigs.remove_contigs()
        self.assertTrue(nx.is_isomorphic(my_graph, remove_contigs.graph))
    
    def test_remove_one_contig(self):
        my_graph = nx.Graph()
        my_graph.add_edges_from([('Contig1', 'Contig2'), ('Contig1', 'Contig3'), ('Contig1', 'Contig4')])
        remove_contigs = RemoveContigs(my_graph)
        expected_graph = nx.Graph()
        expected_graph.add_nodes_from(['Contig1', 'Contig2', 'Contig3'])
        remove_contigs.remove_contigs()
        self.assertTrue(nx.is_isomorphic(expected_graph, remove_contigs.graph))
        