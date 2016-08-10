import unittest
import networkx as nx
from madansi.RefineContigNeighbours import RefineContigNeighbours

class TestRefineContigNeighbours(unittest.TestCase):
    
    def test_count_contig_appearances_empty_list_and_dictionary(self):
        """Tests when both the gene list and the gene to contig dictionary is empty"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance([]), {})
    
    def test_count_contig_appearance_dictionary_empty(self):
        """Tests when there is a gene in the gene list not given in the dictionary"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance(['gene1']), {})
    
    def test_count_contig_appearance_gene_list_empty(self):
        """Tests when the dictionary is not empty but the gene list is"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {'gene1':'Contig1'}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance([]), {})
        
    def test_count_contig_appearance_one_entry(self):
        """Tests when the gene list has one entry present in the dictionary"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2'}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance(['gene1']), {'Contig1':1})
    
    def test_count_contig_appearance_three_entries(self):
        """Tests when there are multiple genes in the same contig present in the list"""
        refine_contig_neighbours = RefineContigNeighbours('','','madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2'}
        self.assertEqual(refine_contig_neighbours.count_contig_appearance(['gene1', 'gene2', 'gene3']), {'Contig1':2, 'Contig2':1})
        
    def test_keep_all_connections_empty_list(self):
        """Tests an empty graph and initial list of neighbouring contigs"""
        filtered_graph = nx.Graph()
        neighbouring_contigs = []
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {}
        refine_contig_neighbours.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours.refined_neighbouring_contigs), sorted(neighbouring_contigs))
    
    def test_keep_one_connection(self):
        """Tests that one connection in neighbouring contigs is kept under refinement"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1', 'gene2'), ('gene2', 'gene5'), ('gene5', 'gene6'), ('gene6', 'gene3'), ('gene3', 'gene4')])
        neighbouring_contigs = [[('Contig1', 'Contig2'), 2, ['gene5', 'gene6']]]
        
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2', 'gene4':'Contig2'}
        refine_contig_neighbours.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours.refined_neighbouring_contigs), sorted(neighbouring_contigs))
    
    def test_keep_three_connections(self):
        """Tests that three connections are preserved under refinement"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'), ('gene4', 'gene5'), ('gene5', 'gene6'), \
                                        ('gene6', 'gene7'), ('gene7', 'gene8'), ('gene8', 'gene9'), ('gene9', 'gene10'), ('gene10', 'gene11')])
        neighbouring_contigs = [[('Contig1', 'Contig2'),1, ['gene3', 'gene4'], [('Contig2', 'Contig3'),2,['gene6','gene7']], \
                                [('Contig3', 'Contig4'),1,['gene9', 'gene10']]]]
        
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = {  'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1', 'gene4':'Contig2', 'gene5':'Contig2', \
                                            'gene8':'Contig3', 'gene9':'Contig3', 'gene10':'Contig4', 'gene11':'Contig4'}
        refine_contig_neighbours.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours.refined_neighbouring_contigs), sorted(neighbouring_contigs))                        
                                                                   
    def test_three_contigs_together(self):
        """Tests when there are three contigs within a small distance of the intersection points"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1', 'gene2'), ('gene3', 'gene4'), ('gene5', 'gene6'), ('gene7','gene2'), ('gene7','gene3'), ('gene7', 'gene5')])
        neighbouring_contigs = [[('Contig1', 'Contig2'),1,['gene7']], [('Contig2', 'Contig3'),1,['gene7']], [('Contig1', 'Contig3'),1,['gene7']]]
        
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = { 'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2', 'gene4':'Contig2', 'gene5':'Contig3', 'gene6':'Contig3'}
        refine_contig_neighbours.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours.refined_neighbouring_contigs), [])
        
    def test_single_gene_from_a_contig(self):
        """Tests when there are two contigs present, one with a single gene on it"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1', 'gene2'), ('gene2', 'gene3')])
        neighbouring_contigs = [[('Contig1', 'Contig2'),1,['gene2', 'gene3']]]
        
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = { 'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2'}
        refine_contig_neighbours.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours.refined_neighbouring_contigs), [])
    
    def test_contig_joins_in_middle(self):
        """Tests when an intersection is found in the case where one end of the contig is closest to the middle of a second"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'), ('gene4', 'gene5'), ('gene2', 'gene6'), ('gene6', 'gene7')])
        neighbouring_contigs = [[('Contig1', 'Contig2'),1,['gene2', 'gene6']]]
        
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/empty_file')
        refine_contig_neighbours.genes = { 'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1', 'gene4':'Contig1', 'gene5':'Contig1', 'gene6':'Contig2', 'gene7':'Contig2'}
        refine_contig_neighbours.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours.refined_neighbouring_contigs), [])
    