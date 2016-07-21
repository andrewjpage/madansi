import unittest
from madansi.DepthFirstSearch import DepthFirstSearch

class TestDepthFirstSearch(unittest.TestCase):
    
    def test__init__(self):
        dfs = DepthFirstSearch('ab','cd')
        self.assertTrue(dfs.graphfile)
        self.assertTrue(dfs.filteredfile)
        
    def test_find_contig(self):
        """Tests that the contig that the gene is in is correctly identified"""
        dfs = DepthFirstSearch('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/gene_present_unittest')
        contig = dfs.find_contig('Contig3')
        self.assertEqual(contig, '7.23.B265.9.cap3_contig')
    
    def test_graph_3_nodes(self):
        """Tests a graph that has 3 nodes which are all marked as present in the lookup table- genes in a chain"""
        dfs = DepthFirstSearch('madansi/tests/data/graph_3_nodes.dot', 'madansi/tests/data/gene_present_unittest')
        list_dfs=[tuple(sorted(tup)) for tup in dfs.depth_first_search()]
        expected_list=[('Contig1','Contig2'),('Contig2','Contig3')]
        self.assertCountEqual(list_dfs, expected_list)
        
    def test_graph_4_nodes(self):
        """Tests a graph that has 4 nodes with one gene not present in the lookup table- genes in a chain"""
        dfs = DepthFirstSearch('madansi/tests/data/graph_4_nodes.dot', 'madansi/tests/data/gene_present_unittest')
        list_dfs=[tuple(sorted(tup))for tup in dfs.depth_first_search()]
        expected_list=[('Contig1','Contig2'),('Contig2','Contig3')]
        self.assertCountEqual(list_dfs, expected_list)
          
        dfs = DepthFirstSearch('madansi/tests/data/graph_4_nodes_2.dot', 'madansi/tests/data/gene_present_unittest')
        list_dfs=[tuple(sorted(tup))for tup in dfs.depth_first_search()]
        expected_list=[('Contig2','Contig3')]
        self.assertCountEqual(list_dfs, expected_list)
        
    def test_graph_5_nodes(self):
        dfs = DepthFirstSearch('madansi/tests/data/graph_5_nodes.dot', 'madansi/tests/data/filtered_data_5_contigs')
        list_dfs=[tuple(sorted(tup)) for tup in dfs.depth_first_search()]
        expected_list=[('Contig1','Contig2'),('Contig4','Contig5')]
        self.assertCountEqual(list_dfs, expected_list)
        
        dfs = DepthFirstSearch('madansi/tests/data/graph_5_nodes_2.dot', 'madansi/tests/data/filtered_data_5_contigs')
        list_dfs=[tuple(sorted(tup)) for tup in dfs.depth_first_search()]
        expected_list=[('Contig1','Contig2'),('Contig2','Contig4')]
        self.assertCountEqual(list_dfs, expected_list)
    
# Should maybe also test add_node_attribute to check that the correct data has been added
        
        