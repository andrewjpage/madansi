import unittest
from madansi.UnusedContigs import UnusedContigs
from madansi.GeneDetector import GeneDetector
import filecmp
import os

class TestUnusedContigs(unittest.TestCase):
    def test_correct_output(self):
        """Tests that the format that the unused contigs are in is correct"""
        output_file = 'output'
        unused_contigs = UnusedContigs('', '', output_file) 
        unused_contigs.unused_contigs_list = ['Contig1', 'Contig2', 'Contig3']
        unused_contigs.output_unused_contigs()
        
        self.assertTrue(filecmp.cmp(output_file, 'madansi/tests/data/unused_contigs', shallow= False))
        os.unlink(output_file)
    
    def test_no_difference(self):
        """Tests the case when all the contigs in the assembly file are also given in the filtered blast hits file"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits')
        output_file = 'output'
        unused_contigs = UnusedContigs(gene_detector, gene_detector.assembly.sequence_names() , output_file)
        unused_contigs.contigs_not_in_filtered_file()
        unused_contigs.output_unused_contigs()
        self.assertTrue(filecmp.cmp(output_file, 'madansi/tests/data/empty_file'))
        os.unlink(output_file)
    
    def test_missing_sequence(self):
        """Tests when there is a difference between the contigs present in the assembly file and those where a gene is present"""
        gene_detector = GeneDetector('madansi/tests/data/assembly.fa', 'madansi/tests/data/test_blast_hits_2')
        output_file = 'output'
        unused_contigs = UnusedContigs(gene_detector, gene_detector.assembly.sequence_names() , output_file)
        unused_contigs.contigs_not_in_filtered_file()
        unused_contigs.output_unused_contigs()
        self.assertTrue(filecmp.cmp(output_file, 'madansi/tests/data/contig3'))
        os.unlink(output_file)
        
            
    