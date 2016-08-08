import unittest
from madansi.GraphParser import GraphParser
import networkx as nx

class TestGraphParse(unittest.TestCase):
    
    def test_initialisation(self):
        graph_parser = GraphParser('madansi/tests/data/test_graph.dot')
        self.assertTrue(graph_parser.graph_file)
        self.assertTrue(graph_parser.graph)
        
    def test_filter_one_heavier(self):
        graph_parser = GraphParser('madansi/tests/data/graph_parser_one_heavier.dot')
        filtered_graph = graph_parser.filter_graph(graph_parser.open_graph_file())
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_parser_expected_one_heavier.dot'))
        self.assertTrue(nx.is_isomorphic(filtered_graph, expected_graph))
    
    def test_filter_no_change(self):
        graph_parser = GraphParser('madansi/tests/data/graph_parser_weighted.dot')
        filtered_graph = graph_parser.filter_graph(graph_parser.open_graph_file())
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_parser_weighted.dot'))
        self.assertTrue(nx.is_isomorphic(filtered_graph, expected_graph))
    
    def test_filter_three_heavier(self):
        graph_parser = GraphParser('madansi/tests/data/graph_parser_three_heavier.dot')
        filtered_graph = graph_parser.filter_graph(graph_parser.open_graph_file())
        expected_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_parser_expected_three_heavier.dot'))
        self.assertTrue(nx.is_isomorphic(filtered_graph, expected_graph))
        