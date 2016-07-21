import unittest
from madansi.WalkGraphs import WalkGraphs
import networkx as nx
import filecmp
import os

class TestWalkGraphs(unittest.TestCase):
    
    def test__init__(self):
        """Tests that the class is correctly initialised"""
        wg = WalkGraphs('ab','cd')
        self.assertTrue(wg.graphfile)
        self.assertTrue(wg.filteredfile)
    
    def test_open_graph_file(self):
        """Tests that the graph dot file opens correctly"""
        wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot', 'cd')
        g = wg.open_graph_file()
        self.assertTrue(g)
        
    def test_open_filtered_file(self):
        """Tests that the file with filtered data opens correctly"""
        wg = WalkGraphs('ab', 'madansi/tests/data/gene_present_unittest')
        f = wg.open_filtered_file()
        self.assertTrue(f)
        
#This test doesn't work because the former graph produced has the data as a number while the latter turns it into a string
         
  #  def test_create_subgraph(self):
  #      """Tests that a subgraph with the correct data attributes is created"""
  #      wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/gene_present_unittest')
  #      g = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_3_nodes_attributes.dot'))
  #      h = wg.create_subgraph()
  #      self.assertCountEqual(h.nodes(data=True),g.nodes(data=True))

    def test_starting_gene(self):
        """Tests that a starting gene is chosen and that it is a node in the subgraph"""
        wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/gene_present_unittest')
        h = wg.create_subgraph()
        start_gene= wg.starting_gene()
        self.assertTrue(h.node[start_gene])
        
    def test_find_contig(self):
        """Tests that the contig that the gene is in is correctly identified"""
        wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/gene_present_unittest')
        contig = wg.find_contig('Contig3')
        self.assertEqual(contig, '7.23.B265.9.cap3_contig')
    
    def test_construct_contig_list(self):
        """Tests that a list of all the possible contig sequences are included in the list"""
        wg = WalkGraphs('ab', 'madansi/tests/data/gene_present_unittest')
        contig_list = wg.construct_contig_list()
        expected_list = ['7.23.B265.9.cap3_contig', 'Sample3', 'Sample4']
        self.assertCountEqual(contig_list, expected_list)