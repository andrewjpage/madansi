import unittest
from madansi.Gene import Gene
import networkx as nx

class TestGene(unittest.TestCase):
    
    def test__basic_setup(self):
        graph = nx.Graph()
        graph.add_node(1)    
        node = graph.nodes()
        my_gene = Gene(1, 1,10,node ,'Contig1', 350)
        self.assertEqual(my_gene.orientation, 1)
        self.assertEqual(my_gene.start, 1)
        self.assertEqual(my_gene.end,10)
        self.assertEqual(my_gene.node, node)
        self.assertEqual(my_gene.contig, 'Contig1')
        self.assertEqual(my_gene.qry_start, 350)