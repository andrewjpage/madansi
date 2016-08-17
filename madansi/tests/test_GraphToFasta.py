import unittest
from madansi.GraphToFasta import GraphToFasta
import os
import filecmp
import networkx as nx


class TestGraphToFasta(unittest.TestCase):
    
    def test_find_contigs_degree_one(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        graph_to_fasta = GraphToFasta('madansi/tests/data/empty_file.fa', graph, 'output.fa', {})
        self.assertEqual(sorted(['Contig1', 'Contig4']), sorted(graph_to_fasta.contigs_degree_one(['Contig1', 'Contig2', 'Contig3', 'Contig4'])))
    
    def test_determine_orientation(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        contig_ends = {'Contig1':{'Contig2': (1,100)}, 'Contig2':{'Contig1':(29, 102), 'Contig3':(307, 240)}, 'Contig3':{'Contig2':(1000,891), 'Contig4':(201, 401)}, 'Contig4':{'Contig3':(300, 49)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/empty_file.fa', graph, 'output.fa', contig_ends)
        self.assertEqual(graph_to_fasta.determine_orientation(['Contig1', 'Contig2'], 'Contig3'), -1)
        self.assertEqual(graph_to_fasta.determine_orientation(['Contig1'], 'Contig2'), 1)
        self.assertEqual(graph_to_fasta.determine_orientation_start_contig('Contig1'), -1)
        self.assertEqual(graph_to_fasta.determine_orientation_end_contig('Contig4'), -1)
    
    def test_walk_contig_graph(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        contig_ends = {'Contig1':{'Contig2': (1,100)}, 'Contig2':{'Contig1':(29, 102), 'Contig3':(307, 240)}, 'Contig3':{'Contig2':(1000, 891), 'Contig4':(201, 401)}, 'Contig4':{'Contig3':(49, 300)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/empty_file.fa', graph, 'output.fa', contig_ends)
        self.assertEqual(['Contig1', 'Contig2', 'Contig3', 'Contig4'], graph_to_fasta.walk_contig_graph(['Contig1', 'Contig2', 'Contig3', 'Contig4']))
    
    def test_walk_contig_graph_contig_high_degree(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4'), ('Contig4', 'Contig5'), ('Contig4', 'Contig7'), ('Contig5', 'Contig6')])
        contig_ends = { 'Contig1':{'Contig2': (1,100)},\
                        'Contig2':{'Contig1':(29, 102), 'Contig3':(307, 240)},\
                        'Contig3':{'Contig2':(1000, 891), 'Contig4':(201, 401)},\
                        'Contig4':{'Contig3':(49, 300), 'Contig5':(102, 30), 'Contig7':(401,209)},\
                        'Contig5':{'Contig4':(1,302), 'Contig6':(1,100)},\
                        'Contig6':{'Contig5':(8,912)},\
                        'Contig7':{'Contig4':(1,20)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/assembly_7_sequences.fa', graph, 'output.fa', contig_ends)
        expected_graph = nx.Graph()
        expected_graph.add_edges_from([('Contig1', 'Contig2'),('Contig2', 'Contig3'), ('Contig4', 'Contig5')])
        self.assertTrue(nx.is_isomorphic(graph_to_fasta.graph, expected_graph))
        self.assertEqual(graph_to_fasta.walk_contig_graph(['Contig1', 'Contig2', 'Contig3']), ['Contig1', 'Contig2', 'Contig3'])
        graph_to_fasta.create_fasta_file_combined_contigs()
        self.assertTrue(filecmp.cmp('madansi/tests/data/combine_contigs.fa', 'output.fa'))
        os.unlink('output.fa')
    
    def test_create_fasta_file(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig3'), ('Contig3', 'Contig2'), ('Contig2', 'Contig4')])
        contig_ends = {'Contig1':{'Contig3': (1, 100)}, 'Contig3':{'Contig1':(29, 102), 'Contig2':(307, 240)}, 'Contig2':{'Contig3':(1000, 891), 'Contig4':(201, 401)}, 'Contig4':{'Contig2':(49, 300)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/assembly_4_sequences.fa', graph, 'output.fa', contig_ends)
        graph_to_fasta.create_fasta_file()
        self.assertTrue(filecmp.cmp('madansi/tests/data/expected_ordered_contigs.fa', 'output.fa'))
        os.unlink('output.fa')
        
    def test_contig_orientation_f_f(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(490,3)}, 'Contig2': {'Contig1':(1,300)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/assembly_4_sequences.fa', graph, 'output.fa', contig_ends)
        self.assertEqual(graph_to_fasta.determine_orientation_start_contig('Contig1'), 1)
        self.assertEqual(graph_to_fasta.determine_orientation_end_contig('Contig2'), 1)
    
    def test_contig_orientation_f_r(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(490,3)}, 'Contig2': {'Contig1':(300,1)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/assembly_4_sequences.fa', graph, 'output.fa', contig_ends)
        self.assertEqual(graph_to_fasta.determine_orientation_start_contig('Contig1'), 1)
        self.assertEqual(graph_to_fasta.determine_orientation_end_contig('Contig2'), -1)
    
    def test_contig_orientation_r_r(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(3,490)}, 'Contig2': {'Contig1':(300,1)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/assembly_4_sequences.fa', graph, 'output.fa', contig_ends)
        self.assertEqual(graph_to_fasta.determine_orientation_start_contig('Contig1'), -1)
        self.assertEqual(graph_to_fasta.determine_orientation_end_contig('Contig2'), -1)
    
    def test_contig_orientation_r_f(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(3,490)}, 'Contig2': {'Contig1':(1,300)}}
        graph_to_fasta = GraphToFasta('madansi/tests/data/assembly_4_sequences.fa', graph, 'output.fa', contig_ends)
        self.assertEqual(graph_to_fasta.determine_orientation_start_contig('Contig1'), -1)
        self.assertEqual(graph_to_fasta.determine_orientation_end_contig('Contig2'), 1)
