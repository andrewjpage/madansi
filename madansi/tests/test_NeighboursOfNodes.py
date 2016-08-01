import unittest
from madansi.NeighboursOfNodes import NeighboursOfNodes
from madansi.GraphParser import GraphParser


class TestNeighboursOfNodes(unittest.TestCase):
    
    def test_initialise(self):
        neighbours_of_nodes = NeighboursOfNodes('') 
        self.assertEqual(neighbours_of_nodes.seen_nodes, [])
    
    def test_find_neighbours_empty_list(self):
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours([])
        self.assertCountEqual(neighbour_list, [])
    
    def test_find_neighbours_one_node(self):
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours(['gene2'])
        self.assertCountEqual(neighbour_list, ['gene1','gene3'])
    
    def test_find_neighbours_three_nodes(self):
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours(['gene3', 'gene4', 'gene5'])
        self.assertCountEqual(neighbour_list, ['gene1', 'gene2','gene6'])
    
    def test_find_neighbours_all_nodes(self):
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours(['gene1','gene2', 'gene3', 'gene4', 'gene5', 'gene6', 'gene7'])
        self.assertCountEqual(neighbour_list, [])