import unittest
from madansi.ContigSpanningTree import ContigSpanningTree
import networkx as nx
import os

class TestContigSpanningTree(unittest.TestCase):
    
    def test_spanning_tree_obtained(self):
        """Tests that a spanning tree is obtained and there is an output file of the graph"""
        graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/test_graph.dot'))
        output_graph = 'output.dot'
        spanning_tree = ContigSpanningTree(graph, output_graph)
        contig_spanning_tree = spanning_tree.construct_spanning_tree()
        self.assertCountEqual(contig_spanning_tree.nodes(), ['gene1','gene2','gene3','gene4','gene5','gene6', 'gene7'])
        self.assertTrue(os.path.isfile(output_graph))
        os.unlink(output_graph)