import unittest
from madansi.WalkGraphs import WalkGraphs

class TestWalkGraphs(unittest.TestCase):
    
    def test__init__(self):
        wg = WalkGraphs('ab','cd')
        self.assertTrue(wg.graphfile)
        self.assertTrue(wg.filteredfile)
    
    def test_graph_3_nodes(self):
        """Tests a graph that has 3 nodes which are all marked as present in the lookup table"""
        wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot', 'madansi/tests/data/gene_present_unittest')
        #start_node=wg.choose_starting_node()
        list_dfs=[tuple(sorted(tup)) for tup in wg.walk_graph()]
        expected_list=[('Contig1','Contig2'),('Contig2','Contig3')]
        self.assertCountEqual(list_dfs, expected_list)