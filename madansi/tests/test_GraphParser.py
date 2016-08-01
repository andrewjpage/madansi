import unittest
from madansi.GraphParser import GraphParser

class TestGraphParse(unittest.TestCase):
    
    def test_initialisation(self):
        graph_parser = GraphParser('madansi/tests/data/test_graph.dot')
        self.assertTrue(graph_parser.graphfile)
        self.assertTrue(graph_parser.graph)
        