import unittest
from madansi.ContigSearching import ContigSearching
from madansi.GeneDetector import GeneDetector
import networkx as nx

class TestContigSearching(unittest.TestCase):
    """Works when one contig has no genes and when there are three contigs"""
    def test_initialisation(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_3_nodes.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        self.assertDictEqual(contig_searching.genes_in_contig_radius, {})
        self.assertCountEqual(contig_searching.neighbouring_contigs, [])
        self.assertEqual(contig_searching.finished_contigs, set())
       
    def test_set_expansion_one_contig_isolated(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_3_nodes.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        self.assertEqual( contig_searching.set_expansion(['gene1','gene2','gene3'], 'Contig1').finished_contigs, set(['Contig1']))
       
    def test_set_expansion_one_contig_not_isolated(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/test_graph.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        self.assertEqual( contig_searching.set_expansion(['gene1','gene2','gene3'], 'Contig1').finished_contigs, set())
    
    def test_expand_all_contigs_two_neighbouring_contigs(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/four_blast_hits')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_4_nodes.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        contig_searching.expand_all_contigs()
        contig_neighbourhoods = contig_searching.neighbouring_contigs
        self.assertTrue(contig_neighbourhoods == [[sorted(('Contig1', 'Contig2')), 1, sorted(['geneA', 'gene1'])]])
            
    def test_expand_all_contigs_three_contigs(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/test_graph.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        contig_searching.expand_all_contigs()
        contig_neighbourhoods = contig_searching.neighbouring_contigs
        self.assertTrue(sorted(contig_neighbourhoods) == sorted([   [sorted(('Contig1', 'Contig2')),1,sorted(['gene1', 'gene4', 'gene3', 'gene5'])] ,\
                                                                    [sorted(('Contig1', 'Contig3')),1,sorted(['gene2', 'gene6'])] ]))
                                                                        
    def test_expand_all_contigs_multiple_iterations(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits_2')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/test_graph_2.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        contig_searching.expand_all_contigs()
        contig_neighbourhoods = contig_searching.neighbouring_contigs
        self.assertTrue(contig_neighbourhoods == [[sorted(('Contig1', 'Contig2')), 3, sorted(['geneB', 'geneD'])]])
    
    def test_expand_all_contigs_three_contigs_multiple_iterations(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits_three_contigs')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/three_contigs_separated.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        contig_searching.expand_all_contigs()
        contig_neighbourhoods = contig_searching.neighbouring_contigs
        self.assertTrue(sorted(contig_neighbourhoods) == sorted([[  sorted(('Contig1', 'Contig2')),2, sorted(['geneA', 'geneD'])],\
                                                                [   sorted(('Contig2', 'Contig3')),2, sorted(['geneB', 'geneE'])],\
                                                                [   sorted(('Contig3', 'Contig1')),2, sorted(['geneC', 'geneF'])]]))
        
    def test_one_contig_dummy_genes(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/one_blast_hit')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/one_contig_dummy_genes.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        contig_searching.expand_all_contigs()
        contig_neighbourhoods = contig_searching.neighbouring_contigs
        self.assertEqual(contig_neighbourhoods, [])
    
    def test_two_separated_sections(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/seven_blast_hits')
        filtered_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/two_separated_sections.dot'))
        contig_searching = ContigSearching(gene_detector, filtered_graph)
        contig_searching.expand_all_contigs()
        contig_neighbourhoods = contig_searching.neighbouring_contigs
        possible_expected_lists = \
        [sorted([('Contig1', 'Contig2', 1), ('Contig3', 'Contig4', 2)]),\
         sorted([('Contig2', 'Contig1', 1), ('Contig3', 'Contig4', 2)]),\
         sorted([('Contig1', 'Contig2', 1), ('Contig4', 'Contig3', 2)]),\
         sorted([('Contig2', 'Contig1', 1), ('Contig4', 'Contig3', 2)])]
        self.assertTrue(sorted(contig_neighbourhoods) == sorted([[  sorted(('Contig1', 'Contig2')), 1, sorted(['gene2', 'gene3'])],\
                                                                 [  sorted(('Contig3', 'Contig4')), 2, sorted(['geneA', 'geneB'])]]))
        

        