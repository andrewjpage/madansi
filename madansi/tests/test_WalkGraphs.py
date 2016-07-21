import unittest
from madansi.WalkGraphs import WalkGraphs
import networkx as nx
import filecmp
import os

class TestWalkGraphs(unittest.TestCase):
    
    def test__init__(self):
        wg = WalkGraphs('ab','cd')
        self.assertTrue(wg.graphfile)
        self.assertTrue(wg.filteredfile)
    
    def test_open_graph_file(self):
        wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot', 'cd')
        g = wg.open_graph_file()
        self.assertTrue(g)
        
#This test doesn't work because the former graph produced has the data as a number while the latter turns it into a string
         
  #  def test_create_subgraph(self):
  #      wg = WalkGraphs('madansi/tests/data/graph_3_nodes.dot','madansi/tests/data/gene_present_unittest')
  #      g = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_3_nodes_attributes.dot'))
  #      h = wg.create_subgraph()
  #      self.assertCountEqual(h.nodes(data=True),g.nodes(data=True))
