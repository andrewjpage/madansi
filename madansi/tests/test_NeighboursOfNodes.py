import unittest
from madansi.NeighboursOfNodes import NeighboursOfNodes
from madansi.GraphParser import GraphParser


class TestNeighboursOfNodes(unittest.TestCase):
    
    def test_initialise(self):
        """Tests initialisation of empty list when given no list"""
        neighbours_of_nodes = NeighboursOfNodes('') 
        self.assertEqual(neighbours_of_nodes.seen_nodes, [])
    
    def test_find_neighbours_empty_list(self):
        """Tests that the method works when given no nodes"""
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours([])
        self.assertCountEqual(neighbour_list, [])
    
    def test_find_neighbours_one_node(self):
        """Tests expected output with one node"""
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours(['gene2'])
        self.assertCountEqual(neighbour_list, ['gene1','gene3', 'gene6'])
    
    def test_find_neighbours_three_nodes(self):
        """Tests that all the neighbours are found when three nodes are included"""
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours(['gene3', 'gene4', 'gene5'])
        self.assertCountEqual(neighbour_list, ['gene1', 'gene2'])
    
    def test_find_neighbours_all_nodes(self):
        """Tests that empty list is returned when all the nodes in the graph are in the list"""
        graph = GraphParser('madansi/tests/data/test_graph.dot').graph
        neighbours_of_nodes = NeighboursOfNodes(graph)
        neighbour_list = neighbours_of_nodes.find_neighbours(['gene1','gene2', 'gene3', 'gene4', 'gene5', 'gene6', 'gene7'])
        self.assertCountEqual(neighbour_list, [])