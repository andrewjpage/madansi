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
        
    def test_construct_contig_list(self):
        """Tests that a list of all the possible contig sequences are included in the list"""
        wg = WalkGraphs('ab', 'madansi/tests/data/gene_present_unittest')
        contig_list = wg.construct_contig_list()
        expected_list = ['7.23.B265.9.cap3_contig', 'Sample3', 'Sample4']
        self.assertCountEqual(contig_list, expected_list)
    
  #  def test_find_neighbors_on_contig(self):
  #      """Tests that the correct neighbors are given on the contig"""
  #      wg = WalkGraphs('madansi/tests/data/graph_4_nodes.dot', 'madansi/tests/data/filtered_data_4_contigs')
  #      contig_neighbor_list = wg.find_neighbors_on_contig('Contig2')
  #      expected_list = ['Contig1','Contig3']
  #      self.assertCountEqual(expected_list, contig_neighbor_list)
    
  #  def test_start_and_end_contigs(self):
  #      """Tests that the correct start and end to the contigs is given"""
  #      wg = WalkGraphs('madansi/tests/data/graph_5_nodes_2','madansi/tests/data/filtered_data_5_contigs')
  #      ends_of_contig = wg.start_and_end_contigs('Contig1')
  #      expected_ends = ['Contig1', 'Contig5']
  #      self.assertCountEqual(ends_of_contig, expected_ends)
  #      
  #      ends_of_contig = wg.start_and_end_contigs('Contig3')
 
    def test_find_ends_of_contig(self):
        """Tests that the ends of the contig have been found correctly"""
        wg = WalkGraphs('madansi/tests/data/graph_5_nodes_2.dot', 'madansi/tests/data/filtered_data_5_contigs_2')
        end_list = wg.find_ends_of_contig('Contig1')
        expected_list = ['Contig1', 'Contig5']
        self.assertCountEqual(end_list,expected_list)
        
        end_list = wg.find_ends_of_contig('Contig3')
        self.assertCountEqual(end_list,expected_list)
        
    def test_order_contigs(self):
        wg = WalkGraphs('madansi/tests/data/graph_order_contigs_1','madansi/tests/data/filtered_data_test_order_contigs')
        output_list = wg.order_contigs()
                