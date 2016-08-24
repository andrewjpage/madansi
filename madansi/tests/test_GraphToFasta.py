import unittest
from madansi.GraphToFasta import GraphToFasta
from madansi.Assembly import Assembly
import os
import filecmp
import networkx as nx


class TestGraphToFasta(unittest.TestCase):
    
    def test_find_contigs_degree_one(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        assembly = Assembly('madansi/tests/data/empty_file.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', {})
        self.assertEqual(sorted(['Contig1', 'Contig4']), sorted(graph_to_fasta.contigs_degree_one(['Contig1', 'Contig2', 'Contig3', 'Contig4'])))
    
    def test_walk_contig_graph(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        contig_ends = {'Contig1':{'Contig2': (1,100)}, 'Contig2':{'Contig1':(29, 102), 'Contig3':(307, 240)}, 'Contig3':{'Contig2':(1000, 891), 'Contig4':(201, 401)}, 'Contig4':{'Contig3':(49, 300)}}
        assembly = Assembly('madansi/tests/data/assembly_4_sequences.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', contig_ends)
        self.assertEqual(['Contig1', 'Contig2', 'Contig3', 'Contig4'], graph_to_fasta.walk_contig_graph(['Contig1', 'Contig2', 'Contig3', 'Contig4']))
    
    def test_walk_contig_graph_contig_high_degree(self):
        contig_ends = { 'Contig1':{'Contig2':(1,100)},\
                        'Contig2':{'Contig1':(29, 102), 'Contig3':(307, 240)},\
                        'Contig3':{'Contig2':(1000, 891), 'Contig4':(201, 401)},\
                        'Contig4':{'Contig3':(49, 300), 'Contig5':(102, 30), 'Contig7':(401,209)},\
                        'Contig5':{'Contig4':(1,302), 'Contig6':(1,100)},\
                        'Contig6':{'Contig5':(8,912)},\
                        'Contig7':{'Contig4':(1,20)}}
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'),('Contig2', 'Contig3'), ('Contig4', 'Contig5')])
        assembly = Assembly('madansi/tests/data/assembly_7_sequences.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', contig_ends)
        
        self.assertEqual(graph_to_fasta.walk_contig_graph(['Contig1', 'Contig2', 'Contig3']), ['Contig1', 'Contig2', 'Contig3'])
        graph_to_fasta.create_fasta_file_combined_contigs()
        self.assertTrue(filecmp.cmp('madansi/tests/data/combine_contigs.fa', 'output.fa'))
        os.unlink('output.fa')
    
    def test_3_cycle_f_f_f(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig1')])
        contig_ends = {'Contig1':{'Contig2':(100,1), 'Contig3':(1,100)}, 'Contig2':{'Contig1':(1,100), 'Contig3':(100,1)}, 'Contig3':{'Contig1':(100,1), 'Contig2':(1,100)}}
        assembly = Assembly('madansi/tests/data/assembly_7_sequences.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', contig_ends)
        
        
        self.assertEqual(graph_to_fasta.walk_contig_graph(['Contig2', 'Contig1', 'Contig3']), ['Contig1', 'Contig2', 'Contig3'])
        self.assertDictEqual(graph_to_fasta.contig_orientation, {'Contig1':1, 'Contig2':1, 'Contig3':1})
        
    def test_3_cycle_r_r_r(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig1')])
        contig_ends = {'Contig1':{'Contig2':(1,100), 'Contig3':(100,1)}, 'Contig2':{'Contig1':(100,1), 'Contig3':(1,100)}, 'Contig3':{'Contig1':(1,100), 'Contig2':(100,1)}}
        assembly = Assembly('madansi/tests/data/assembly_7_sequences.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', contig_ends)
        
        self.assertEqual(graph_to_fasta.walk_contig_graph(['Contig2', 'Contig1', 'Contig3']), ['Contig1', 'Contig3', 'Contig2'])
        self.assertDictEqual(graph_to_fasta.contig_orientation, {'Contig1':1, 'Contig2':1, 'Contig3':1})
        
    def test_3_cycle_f_r_f(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig1')])
        contig_ends = {'Contig1':{'Contig2':(100,1), 'Contig3':(1,100)}, 'Contig2':{'Contig1':(100,1), 'Contig3':(1,100)}, 'Contig3':{'Contig1':(100,1), 'Contig2':(1,100)}}
        assembly = Assembly('madansi/tests/data/assembly_7_sequences.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', contig_ends)
        
        self.assertEqual(graph_to_fasta.walk_contig_graph(['Contig2', 'Contig1', 'Contig3']), ['Contig1', 'Contig2', 'Contig3'])
        self.assertDictEqual(graph_to_fasta.contig_orientation, {'Contig1':1, 'Contig2':-1, 'Contig3':1})
        
    def test_4_cycle(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4'), ('Contig4', 'Contig1')])
        contig_ends = {'Contig1':{'Contig2':(100,1), 'Contig4':(1,100)}, 'Contig2':{'Contig1':(1,100), 'Contig3':(100,1)}, 'Contig3':{'Contig4':(100,1), 'Contig2':(1,100)}, 'Contig4':{'Contig1':(100,1), 'Contig3':(1,100)}}
        assembly = Assembly('madansi/tests/data/assembly_7_sequences.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', contig_ends)
        
        
        self.assertEqual(graph_to_fasta.walk_contig_graph(['Contig2', 'Contig1', 'Contig4','Contig3']), ['Contig1', 'Contig2', 'Contig3', 'Contig4'])
        self.assertDictEqual(graph_to_fasta.contig_orientation, {'Contig1':1, 'Contig2':1, 'Contig3':1, 'Contig4':1})
    
    def test_multiple_components(self):
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig4', 'Contig5')])
        contig_ends = {'Contig1':{'Contig2':(100,1)}, 'Contig2':{'Contig1':(1, 100), 'Contig3':(100,1)}, 'Contig3':{'Contig2':(1,100)},'Contig4':{'Contig5':(100,1)}, 'Contig5':{'Contig4':(1,100)}}
        
        assembly = Assembly('madansi/tests/data/assembly_7_sequences.fa')
        assembly.sequence_names()
        graph_to_fasta = GraphToFasta(assembly.sequences, graph, 'output.fa', contig_ends)
        
        graph_to_fasta.create_fasta_file_combined_contigs()
        self.assertTrue(filecmp.cmp('madansi/tests/data/combine_contigs_multiple_components.fa', 'output.fa'))
        os.unlink('output.fa')
        