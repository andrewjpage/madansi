import unittest
import networkx as nx
from madansi.RefineContigNeighbours import RefineContigNeighbours
from madansi.GeneDetector import GeneDetector

class TestRefineContigNeighboursOrientation(unittest.TestCase):
    
    def test_finds_orientation_two_contigs(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'geneA'),\
                                        ('geneA', 'geneB'), ('geneB', 'gene4'), ('gene4', 'gene5')])
        neighbouring_contigs = [[('Contig1', 'Contig2'), 2, ['geneA', 'geneB']]]
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file')
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file', gene_detector)
        
        refine_contig_neighbours_object.contigs = { 'Contig1':{'gene1':None, 'gene2':None, 'gene3':None},\
                                                    'Contig2':{'gene4':None, 'gene5':None}}
                                                    
        refine_contig_neighbours_object.genes = {   'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1',\
                                                    'gene4':'Contig2', 'gene5':'Contig2'}
                                                
        expected_dict = {'Contig1': {'Contig2':(980,490)}, 'Contig2': {'Contig1':(2,301)}}
        refine_contig_neighbours_object.ends_of_contigs()
        self.assertDictEqual(refine_contig_neighbours_object.contig_ends, expected_dict) 
    
    def test_add_to_contig_appearances(self):
        filtered_graph = nx.Graph()
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file')
        neighbouring_contigs = []
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file', gene_detector)
        self.assertDictEqual(refine_contig_neighbours_object.add_to_contig_appearance('gene1', {}), {'Contig1':[1,['gene1',None]]})
        self.assertDictEqual(refine_contig_neighbours_object.add_to_contig_appearance('gene1', {'Contig1':[0, [None, None]]}), {'Contig1':[1,['gene1',None]]})
        self.assertDictEqual(refine_contig_neighbours_object.add_to_contig_appearance('gene1', {'Contig1':[1, ['gene2', None]]}), {'Contig1':[2,['gene2', 'gene1']]})
        self.assertDictEqual(refine_contig_neighbours_object.add_to_contig_appearance('gene1', {'Contig1':[2, ['gene2', 'gene3']]}), {'Contig1':[3,['gene2', 'gene3']]})
        
    def test_orientation_of_contigs(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'),\
                                        ('gene4', 'gene5'), ('gene5', 'gene6'), ('gene6', 'gene7'),\
                                        ('gene7', 'gene8')])
                                
        neighbouring_contigs = [[('Contig1', 'Contig2'),1, ['gene3','gene4']], [('Contig2', 'Contig3'),1,['gene5','gene6']]]
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file' )
        
        refine_contig_neighbours_object =  RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file', gene_detector)
        refine_contig_neighbours_object.genes = {   'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1',\
                                    'gene4':'Contig2', 'gene5':'Contig2', 'gene6':'Contig3', \
                                    'gene7':'Contig3', 'gene8':'Contig3'}
        refine_contig_neighbours_object.ends_of_contigs()
        self.assertDictEqual(refine_contig_neighbours_object.contig_ends, {'Contig1': {'Contig2':(980,490)}, 'Contig2':{'Contig1':(2,301), 'Contig3':(301,2)}, 'Contig3':{'Contig2':(1,250)}})
    
    def test_orientation_further_separation(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'geneA'),\
                                        ('geneA', 'geneB'), ('geneB', 'gene4'), ('gene4', 'gene5'),\
                                        ('gene5', 'geneC'), ('geneC', 'gene6'), ('gene6', 'gene7'),\
                                        ('gene7', 'gene8')])
        neighbouring_contigs = [[('Contig1', 'Contig2'),1, ['geneA','geneB']], [('Contig2', 'Contig3'),1,['geneC']]]
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file' )
        
        refine_contig_neighbours_object =  RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/refine_contig_neighbours_8_blast_hits_file', gene_detector)
        refine_contig_neighbours_object.genes = {   'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1',\
                                    'gene4':'Contig2', 'gene5':'Contig2', 'gene6':'Contig3', \
                                    'gene7':'Contig3', 'gene8':'Contig3', 'geneA':None, 'geneB':None, 'geneC':None}
        refine_contig_neighbours_object.ends_of_contigs()
        self.assertDictEqual(refine_contig_neighbours_object.contig_ends, {'Contig1': {'Contig2':(980,490)}, 'Contig2':{'Contig1':(2,301), 'Contig3':(301,2)}, 'Contig3':{'Contig2':(1,250)}})
                                                
                                                
                                                
                                                
                                                
                                                
                                                
    