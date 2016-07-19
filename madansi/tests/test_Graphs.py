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
        g=Graphs('madansi/tests/data/graph_3_nodes.dot')
        g.open_graph()
        self.assertTrue(g)
        
        #Also test if given the wrong format for a graph or if the graph file does not exist
        
    def test_comparison_graph_lookuptable(self):
        pass