import unittest
from madansi.ContigSearching import ContigSearching
from madansi.GeneDetector import GeneDetector
import networkx as nx

class TestContigSearching(unittest.TestCase):
    
    def test_initialisation(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        filtered_graph = 'madansi/tests/data/test_graph.dot'
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        self.assertDictEqual(contig_searching.genes_in_contig_radius, {})
        self.assertCountEqual(contig_searching.neighbouring_contigs, [])
        self.assertEqual(contig_searching.finished_contigs, set())
        
    def test_set_expansion_one_contig(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_3_nodes.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        self.assertEqual(contig_searching.set_expansion(['gene1','gene2','gene3'], 'Contig1').finished_contigs, set('Contig1'))