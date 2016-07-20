import unittest
from madansi.Graphs import Graphs, Error

class TestGraphs(unittest.TestCase):
    
    def test__init__(self):
        """Tests whether input string is recognised"""
        g = Graphs('ab','cd')
        self.assertTrue(g.graphfile)
        self.assertTrue(g.filteredfile)
    
    def test_open_graph(self):
        """Tests that the dot file can be opened for reading"""
        g=Graphs('madansi/tests/data/graph_3_nodes.dot','cd')
        g.open_graph()
        self.assertTrue(g)
        
    def test_same_genes(self):
        """Tests graph that shows exactly the same genes present as the dictionary"""
        g=Graphs('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/gene_present_unittest')
        obtained_list=g.extra_genes_in_graph()
        expected_list = []
        self.assertCountEqual(obtained_list, expected_list)
        obtained_list=g.genes_not_in_graph()
        expected_list = []
        self.assertCountEqual(obtained_list, expected_list)
                
    def test_extra_genes_in_graph(self):
        """Tests graph that has additional genes to those given in the dictionary"""
        g=Graphs('madansi/tests/data/graph_4_nodes.dot','madansi/tests/data/gene_present_unittest')
        obtained_list=g.extra_genes_in_graph()
        expected_list = ['Contig4']
        self.assertCountEqual(obtained_list, expected_list)
    
    def test_genes_not_in_graph(self):
        """Tests graph that does not have all the genes listed in the dictionary"""
        g=Graphs('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/graphs_unittest_4_genes_present')
        obtained_list=g.genes_not_in_graph()
        expected_list=['Contig4']
        self.assertCountEqual(obtained_list, expected_list)
        

        
        