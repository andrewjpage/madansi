import networkx as nx
import unittest
from madansi.GeneDetector import GeneDetector
from madansi.ProduceOrderedContigGraph import ProduceOrderedContigGraph

class TestProduceOrderedContigGraph(unittest.TestCase):
    
    def test_two_contigs_unconnected(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1', 'gene2'), ('gene2', 'gene3'),('gene4', 'gene5')])              
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file' )
        filtered_blast_hits_file = 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file'
        
        expected_ordered_contig_graph = nx.Graph()
        
        produced_ordered_contig_graph = ProduceOrderedContigGraph(gene_detector, filtered_graph, filtered_blast_hits_file,  {'Contig1':['','',1000], 'Contig2':['', '', 800], 'Contig3':['','',2000]})
        self.assertTrue(nx.is_isomorphic(expected_ordered_contig_graph, produced_ordered_contig_graph.produce_ordered_contig_graph()))
        
    def test_two_contigs_separated_by_one(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'geneA'),\
                                        ('geneA', 'gene6'), ('gene6', 'gene7'), ('gene7', 'gene8')])
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file' )
        filtered_blast_hits_file = 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file'
                                    
        expected_ordered_contig_graph = nx.Graph()
        expected_ordered_contig_graph.add_edge('Contig1', 'Contig3')
        
        produced_ordered_contig_graph = ProduceOrderedContigGraph(gene_detector, filtered_graph, filtered_blast_hits_file,  {'Contig1':['','',1000], 'Contig2':['', '', 800], 'Contig3':['','',2000]})    
        self.assertTrue(nx.is_isomorphic(expected_ordered_contig_graph, produced_ordered_contig_graph.produce_ordered_contig_graph()))
                               
    
    def test_two_contigs_separated_by_two(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'geneA'),\
                                        ('geneA', 'geneB'), ('geneB', 'gene6'), ('gene6', 'gene7')])
                                        
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file' )
        filtered_blast_hits_file = 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file'
        
        expected_ordered_contig_graph = nx.Graph()
        expected_ordered_contig_graph.add_edge('Contig1', 'Contig3')
        
        produced_ordered_contig_graph = ProduceOrderedContigGraph(gene_detector, filtered_graph, filtered_blast_hits_file, {'Contig1':['','',1000], 'Contig2':['', '', 800], 'Contig3':['','',2000]})    
        self.assertTrue(nx.is_isomorphic(expected_ordered_contig_graph, produced_ordered_contig_graph.produce_ordered_contig_graph()))
        
        
    def test_two_contigs_separated_by_three(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'geneA'),\
                                        ('geneA',    'geneB'), ('geneB', 'geneC'), ('geneC', 'gene4'),\
                                        ('gene4', 'gene5')])
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file' )
        filtered_blast_hits_file = 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file'
                                        
        expected_ordered_contig_graph = nx.Graph()
        expected_ordered_contig_graph.add_edge('Contig1', 'Contig2')
                                        
        produced_ordered_contig_graph = ProduceOrderedContigGraph(gene_detector, filtered_graph, filtered_blast_hits_file, {'Contig1':['','',1000], 'Contig2':['', '', 800], 'Contig3':['','',2000]})
        self.assertTrue(nx.is_isomorphic(expected_ordered_contig_graph, produced_ordered_contig_graph.produce_ordered_contig_graph()))
       
    
