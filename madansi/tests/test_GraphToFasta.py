import unittest
from madansi.GraphToFasta import GraphToFasta
import os
import filecmp
import networkx as nx


class TestGraphToFasta(unittest.TestCase):
    
    def test_find_contigs_degree_one(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        graph_to_fasta = GraphToFasta('', graph, '', {})
        self.assertEqual(sorted(['Contig1', 'Contig4']), sorted(graph_to_fasta.find_contigs_degree_one(['Contig1', 'Contig2', 'Contig3', 'Contig4'])))
    
    def test_determine_orientation(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        contig_ends = {'Contig1':{'Contig2': 1}, 'Contig2':{'Contig1':29, 'Contig3':307}, 'Contig3':{'Contig2':1000, 'Contig4':201}, 'Contig4':{'Contig3':49}}
        graph_to_fasta = GraphToFasta('', graph, '', contig_ends)
        self.assertEqual(graph_to_fasta.determine_orientation(['Contig1', 'Contig2'], 'Contig3'), -1)
        self.assertEqual(graph_to_fasta.determine_orientation(['Contig1'], 'Contig2'), 1)
        self.assertEqual(graph_to_fasta.determine_orientation([], 'Contig1'), 1)
        self.assertEqual(graph_to_fasta.determine_orientation(['Contig1', 'Contig2', 'Contig3'], 'Contig4'), 1)
    
    def test_walk_contig_graph(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        contig_ends = {'Contig1':{'Contig2': 1}, 'Contig2':{'Contig1':29, 'Contig3':307}, 'Contig3':{'Contig2':1000, 'Contig4':201}, 'Contig4':{'Contig3':49}}
        graph_to_fasta = GraphToFasta('', graph, '', contig_ends)
        self.assertEqual(['Contig1', 'Contig2', 'Contig3', 'Contig4'], graph_to_fasta.walk_contig_graph(['Contig1', 'Contig2', 'Contig3', 'Contig4']))
    
    def test_create_fasta_file(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig3'), ('Contig3', 'Contig2'), ('Contig2', 'Contig4')])
        contig_ends = {'Contig1':{'Contig3': 1}, 'Contig3':{'Contig1':29, 'Contig2':307}, 'Contig2':{'Contig3':1000, 'Contig4':201}, 'Contig4':{'Contig2':49}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/assembly_4_sequences.fa', graph, 'output.fa', contig_ends)
        graph_to_fasta.create_fasta_file()
        self.assertTrue(filecmp.cmp('madansi/tests/data/expected_ordered_contigs.fa', 'output.fa'))
        
        