import unittest
import networkx as nx
from madansi.ContigOrientation  import ContigOrientation
from madansi.GeneDetector       import GeneDetector

class TestContigOrientation(unittest.TestCase):
    
    def test_find_contig_orientation(self):
        gene_detector               = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/four_blast_hits')
        contig_orientation_object   = ContigOrientation(nx.Graph(), gene_detector)
        
        self.assertEqual(contig_orientation_object.find_contig_orientation('Contig1'), -1)
        self.assertEqual(contig_orientation_object.find_contig_orientation('Contig2'), 1)
    
    def test_orientate_edges_one_edge(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/seven_blast_hits')
        
        contig_graph_1 = nx.Graph()
        contig_graph_2 = nx.Graph()
        contig_graph_3 = nx.Graph()
        contig_graph_4 = nx.Graph()
        contig_graph_1.add_edge('Contig2', 'Contig4', weight = 0.5)
        contig_graph_2.add_edge('Contig2', 'Contig3', weight = 0.5)
        contig_graph_3.add_edge('Contig1', 'Contig2', weight = 0.5)
        contig_graph_4.add_edge('Contig1', 'Contig3', weight = 0.5)
        
        contig_orientation_1_object   = ContigOrientation(contig_graph_1, gene_detector)
        contig_orientation_1_object.repeat_all_connected_components()
        contig_orientation_2_object   = ContigOrientation(contig_graph_2, gene_detector)
        contig_orientation_2_object.repeat_all_connected_components()
        contig_orientation_3_object   = ContigOrientation(contig_graph_3, gene_detector)
        contig_orientation_3_object.repeat_all_connected_components()
        contig_orientation_4_object   = ContigOrientation(contig_graph_4, gene_detector)
        contig_orientation_4_object.repeat_all_connected_components()
        
        self.assertEqual(contig_orientation_1_object.contig_graph.edge['Contig2']['Contig4']['weight'], 1)
        self.assertEqual(contig_orientation_2_object.contig_graph.edge['Contig2']['Contig3']['weight'], 2)
        self.assertEqual(contig_orientation_3_object.contig_graph.edge['Contig1']['Contig2']['weight'], 3)
        self.assertEqual(contig_orientation_4_object.contig_graph.edge['Contig1']['Contig3']['weight'], 4)
    
    def test_orientate_edges_three_edges(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/seven_blast_hits')
        
        contig_graph = nx.Graph()
        contig_graph.add_edge('Contig2', 'Contig4', weight = 0.5)
        contig_graph.add_edge('Contig4', 'Contig1', weight = 0.5)
        contig_graph.add_edge('Contig1', 'Contig3', weight = 0.5)
        
        contig_orientation_object = ContigOrientation(contig_graph, gene_detector)
        contig_orientation_object.repeat_all_connected_components()
        
        self.assertEqual(contig_orientation_object.contig_graph.edge['Contig2']['Contig4']['weight'], 1)
        self.assertEqual(contig_orientation_object.contig_graph.edge['Contig4']['Contig1']['weight'], 2)
        self.assertEqual(contig_orientation_object.contig_graph.edge['Contig1']['Contig3']['weight'], 4)
        