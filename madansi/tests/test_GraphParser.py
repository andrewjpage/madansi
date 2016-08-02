import unittest
from madansi.GraphParser import GraphParser
import networkx as nx

class TestGraphParse(unittest.TestCase):
    
    def test_initialisation(self):
        graph_parser = GraphParser('madansi/tests/data/test_graph.dot')
        self.assertTrue(graph_parser.graph_file)
        self.assertTrue(graph_parser.graph)
        
    def test_filter(self):
        graph_parser = GraphParser('madansi/tests/data/unfiltered_graph.dot')
        filtered_graph = graph_parser.filter_graph(graph_parser.open_graph_file())
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/filtered_graph.dot'))
        self.assertTrue(nx.is_isomorphic(filtered_graph, expected_graph))
    
    def test_filter_no_change(self):
        graph_parser = GraphParser('madansi/tests/data/test_graph.dot')
        filtered_graph = graph_parser.filter_graph(graph_parser.open_graph_file())
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/test_graph.dot'))
        self.assertTrue(nx.is_isomorphic(filtered_graph, expected_graph))