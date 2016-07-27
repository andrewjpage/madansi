import unittest
from madansi.GenerateGraph import GenerateGraph
import networkx as nx
import filecmp
import os

class TestGenerateGraph(unittest.TestCase):
    
    def test__init__(self):
        gg = GenerateGraph('ab','cd','ef','gh')
        self.assertTrue(gg.graphfile)
        self.assertTrue(gg.filteredfile)
        self.assertTrue(gg.outputgraphfile)
        self.assertTrue(gg.unused_sequence_file)
    
    def test_open_graph_file(self):
        gg = GenerateGraph('madansi/tests/data/graph_3_nodes.dot', 'cd', 'ef', 'gh')
        g = gg.open_graph_file()
        self.assertTrue(g)
    
    def tests_create_dictionaries(self):
        gg = GenerateGraph('ab', 'madansi/tests/data/gene_present_unittest', 'ef', 'gh')
        index_file = gg.create_index_file()
        expected_dict = {'gene2':0, 'gene3':1, 'gene1':2, 'gene4':3}
        self.assertDictEqual(expected_dict,index_file)
        
        gene_present = gg.create_gene_present_dict()
        expected_dict = {'gene2': True, 'gene3':True, 'gene1':True, 'gene4':False}
        self.assertDictEqual(gene_present, expected_dict)
        
        gene_sequence = gg.create_gene_sequence_dict()
        expected_dict = {'gene2':'7.23.B265.9.cap3_contig', 'gene3':'7.23.B265.9.cap3_contig', 'gene1':'Sample3', 'gene4':'Sample4'}
        self.assertDictEqual(gene_sequence, expected_dict)
        
        gene_orientation = gg.create_gene_orientation_dict()
        expected_dict = {'gene2': False, 'gene3': True, 'gene1': False, 'gene4':True}
        self.assertDictEqual(gene_orientation, expected_dict)
    