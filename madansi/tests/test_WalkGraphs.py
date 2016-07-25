import unittest
from madansi.WalkGraphs import WalkGraphs
import networkx as nx
import filecmp
import os

class TestWalkGraphs(unittest.TestCase):
    
    def test__init__(self):
        """Tests that the class is correctly initialised"""
        wg = WalkGraphs('ab','cd', 'ef')
        self.assertTrue(wg.graphfile)
        self.assertTrue(wg.filteredfile)
        self.assertTrue(wg.outputgraphfile)
    
    def test_open_graph_file(self):
        """Tests that the graph dot file opens correctly"""
        wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot', 'cd', 'ef')
        g = wg.open_graph_file()
        self.assertTrue(g)
             
    def test_starting_gene(self):
        """Tests that a starting gene is chosen and that it is a node in the subgraph"""
        wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/gene_present_unittest','ef')
        h = wg.create_subgraph()
        start_gene= wg.starting_gene()
        self.assertTrue(h.node[start_gene])
        
    def test_construct_sequence_list(self):
        """Tests that a list of all the possible sequences are included in the list"""
        wg = WalkGraphs('madansi/tests/data/graph_4_nodes.dot', 'madansi/tests/data/gene_present_unittest','ef')
        sequence_list = wg.construct_sequence_list()
        expected_list = ['7.23.B265.9.cap3_contig', 'Sample3', 'Sample4']
        self.assertCountEqual(sequence_list, expected_list)
    
    def test_find_ends_of_sequence(self):
        """Tests that the ends of the sequence have been found correctly"""
        wg = WalkGraphs('madansi/tests/data/graph_5_nodes_2.dot', 'madansi/tests/data/filtered_data_5_contigs_2', 'ef')
        end_list = wg.find_ends_of_sequence('gene1')
        expected_list = ['gene1', 'gene5']
        self.assertCountEqual(end_list,expected_list)
        
        end_list = wg.find_ends_of_sequence('gene3')
        self.assertCountEqual(end_list,expected_list)
        
    def test_closest_gene(self):
        """Tests that the closest gene on a different sequence to the given end gene is found - linear graph"""
        wg = WalkGraphs('madansi/tests/data/graph_order_sequences.dot','madansi/tests/data/filtered_data_test_order_sequences', 'ef')
        
        gene = wg.closest_gene('gene1')
        self.assertEqual(gene, 'gene6')
        
        gene = wg.closest_gene('gene3')
        self.assertEqual(gene, 'gene6')
        
        gene = wg.closest_gene('gene6')
        self.assertEqual(gene, 'gene3')
        
        gene = wg.closest_gene('gene8')
        self.assertEqual(gene, 'gene3')
       
    def test_order_sequences_cycle(self):
        """Tests that the closest gene on a different sequence to the given end gene is found - cyclic graph"""
        wg = WalkGraphs('madansi/tests/data/graph_order_cycle.dot', 'madansi/tests/data/filtered_data_test_order_sequences','ef')
        
        gene = wg.closest_gene('gene6')
        self.assertEqual(gene,'gene3')
          
        gene = wg.closest_gene('gene3')
        self.assertEqual(gene, 'gene6')
    
        gene = wg.closest_gene('gene8')
        self.assertEqual(gene, 'gene1')
  
        gene = wg.closest_gene('gene1')
        self.assertEqual(gene, 'gene8')
        
    def test_dictionary_pairs_closest_genes(self):
        """Tests that the correct dictionary of closest genes is obtained - linear graph"""
        wg = WalkGraphs('madansi/tests/data/graph_order_sequences.dot','madansi/tests/data/filtered_data_test_order_sequences', 'ef')
        end_genes_dict = wg.dictionary_pairs_closest_genes()
        expected_dictionary = {'gene1': None, 'gene3':'gene6', 'gene6':'gene3', 'gene8':None}   
        self.assertEqual(expected_dictionary, end_genes_dict)
   
    def test_dictionary_pairs_closest_genes_cycle(self):
        """Tests that the correct dictionary is obtained - cyclic graph"""
        wg = WalkGraphs('madansi/tests/data/graph_order_cycle.dot', 'madansi/tests/data/filtered_data_test_order_sequences', 'ef')
        end_genes_dict = wg.dictionary_pairs_closest_genes()
        expected_dictionary = {'gene1':'gene8', 'gene3':None, 'gene6':None, 'gene8':'gene1'}
        self.assertEqual(expected_dictionary, end_genes_dict)
                
    def test_ordering_sequences(self):
        
        wg = WalkGraphs('madansi/tests/data/graph_order_sequences.dot', 'madansi/tests/data/filtered_data_test_order_sequences', 'ef')
        order_genes_visited = wg.ordering_sequences()
        expected_order = ['gene1', 'gene3', 'gene6', 'gene8']
        self.assertTrue(order_genes_visited == expected_order or order_genes_visited==expected_order[::-1])
        
        wg = WalkGraphs('madansi/tests/data/graph_order_cycle.dot', 'madansi/tests/data/filtered_data_test_order_sequences', 'ef')
        order_genes_visited = wg.ordering_sequences()
        expected_order = ['gene6', 'gene8', 'gene1',  'gene3']
        self.assertTrue(order_genes_visited == expected_order or order_genes_visited == expected_order[::-1])
    
 
  
    def test_create_linear_subgraph(self):
        wg = WalkGraphs('madansi/tests/data/graph_order_cycle.dot', 'madansi/tests/data/filtered_data_test_order_sequences', 'output.dot')
        wg.create_linear_subgraph()
        os.unlink('output.dot')
        
        #Need to compare the list of unused sequences with expected and that the produced graph is the same as the expected 

        