import unittest
from madansi.UnusedContigs import UnusedContigs
from madansi.GeneDetector import GeneDetector
import filecmp
import networkx as nx
import os

class TestUnusedContigs(unittest.TestCase):
    def test_correct_output(self):
        """Tests that the format that the unused contigs are in is correct"""
        output_file = 'output.fa'
        unused_contigs = UnusedContigs('', output_file, 'madansi/tests/data/assembly.fa') 
        unused_contigs.unused_contigs_list = ['Contig1', 'Contig2', 'Contig3']
        unused_contigs.add_unused_contigs_to_end()
        
        self.assertTrue(filecmp.cmp(output_file, 'madansi/tests/data/assembly.fa', shallow= False))
        os.unlink(output_file)
    
    def test_no_difference(self):
        """Tests the case when all the contigs in the assembly file are also given in the filtered blast hits file"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        output_file = 'output.fa'
        unused_contigs = UnusedContigs(gene_detector, output_file, 'madansi/tests/data/assembly.fa')
        unused_contigs.contigs_not_in_filtered_file()
        unused_contigs.add_unused_contigs_to_end()
        self.assertTrue(filecmp.cmp(output_file, 'madansi/tests/data/empty_file'))
        os.unlink(output_file)
    
    def test_missing_sequence(self):
        """Tests when there is a difference between the contigs present in the assembly file and those where a gene is present"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits_2')
        output_file = 'output'
        unused_contigs = UnusedContigs(gene_detector, output_file, 'madansi/tests/data/assembly.fa')
        unused_contigs.contigs_not_in_filtered_file()
        unused_contigs.add_unused_contigs_to_end()
        self.assertTrue(filecmp.cmp(output_file, 'madansi/tests/data/contig3'))
        os.unlink(output_file)
    
    def test_compare_graphs(self):
        gene_detector = GeneDetector('madansi/tests/data/assembly_7_sequences.fa', 'madansi/tests/data/test_blast_hits_2')
        output_file = 'output'
        unused_contigs = UnusedContigs(gene_detector, output_file, 'madansi/tests/data/assembly_7_sequences.fa')
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        filtered_graph.add_node('Contig7')
        unused_contigs.contigs_not_in_filtered_graph(filtered_graph)
        unused_contigs.add_unused_contigs_to_end()
        self.assertTrue(filecmp.cmp(output_file, 'madansi/tests/data/test_unused_contigs_in_graph'))
        os.unlink(output_file)
            
    