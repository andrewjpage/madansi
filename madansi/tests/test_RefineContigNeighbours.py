import unittest
import networkx as nx
from madansi.RefineContigNeighbours import RefineContigNeighbours
from madansi.GeneDetector import GeneDetector

class TestRefineContigNeighbours(unittest.TestCase):
    
    def test_keep_all_connections_empty_list(self):
        """Tests an empty graph and initial list of neighbouring contigs"""
        filtered_graph = nx.Graph()
        neighbouring_contigs = []
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/empty_file' )
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/empty_file',gene_detector, {})
        refine_contig_neighbours_object.genes = {}
        refine_contig_neighbours_object.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours_object.refined_neighbouring_contigs), sorted(neighbouring_contigs))
    
    def test_keep_one_connection(self):
        """Tests that one connection in neighbouring contigs is kept under refinement"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene5'), ('gene5', 'gene6'),\
                                        ('gene6', 'gene3'), ('gene3', 'gene4')])
        neighbouring_contigs = [[('Contig1', 'Contig2'), 2, ['gene5', 'gene6']]]
        
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/empty_file' )
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/empty_file',gene_detector, {})
        refine_contig_neighbours_object.genes = {'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2', 'gene4':'Contig2'}
        refine_contig_neighbours_object.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours_object.refined_neighbouring_contigs), sorted(neighbouring_contigs))
    
    def test_keep_two_connections(self):
       """Tests that two connections are preserved under refinement"""
       filtered_graph = nx.Graph()
       filtered_graph.add_edges_from([  ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'),\
                                        ('gene4', 'gene5'), ('gene5', 'gene6'), ('gene6', 'gene7'),\
                                        ('gene7', 'gene8'), ('gene8', 'gene9'), ('gene9', 'gene10'),\
                                        ('gene10', 'gene11'), ('gene11', 'gene12'), ('gene12', 'gene13')])
                                        
       neighbouring_contigs = [[('Contig1', 'Contig2'),1, ['gene4']], [('Contig2', 'Contig3'),2,['gene8','gene9']]]
       gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/empty_file' )
       
       refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs,filtered_graph,'madansi/tests/data/empty_file', gene_detector, {})
       refine_contig_neighbours_object.genes = {   'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1',\
                                            'gene5':'Contig2', 'gene6':'Contig2', 'gene7':'Contig2', \
                                            'gene10':'Contig3', 'gene11':'Contig3', 'gene11':'Contig3',\
                                            'gene12':'Contig3'}
       refine_contig_neighbours_object.refine_contig_neighbours()
       self.assertEqual(sorted(refine_contig_neighbours_object.refined_neighbouring_contigs), sorted(neighbouring_contigs))                        
                                                                  
    def test_three_contigs_together(self):
        """Tests when there are three contigs within a small distance of the intersection points"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene3', 'gene4'), ('gene5', 'gene6'),\
                                        ('gene7','gene2'), ('gene7','gene3'), ('gene7', 'gene5')])
                                        
        neighbouring_contigs =  [[('Contig1', 'Contig2'),1,['gene7']], [('Contig2', 'Contig3'),1,['gene7']],\
                                [('Contig1', 'Contig3'),1,['gene7']]]
                                
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/empty_file' )
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/empty_file',gene_detector, {})
        refine_contig_neighbours_object.genes = {  'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2',\
                                            'gene4':'Contig2', 'gene5':'Contig3', 'gene6':'Contig3'}
        refine_contig_neighbours_object.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours_object.refined_neighbouring_contigs), [])
       
    def test_single_gene_from_a_contig(self):
        """Tests when there are two contigs present, one with a single gene on it when other genes from the same contig are elsewhere in the graph"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1', 'gene2'), ('gene2', 'gene3'), ('gene4', 'gene5')])
        neighbouring_contigs = [[('Contig1', 'Contig2'),1,['gene2', 'gene3']]]
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_5_blast_hits_file' )
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/refine_contig_neighbours_5_blast_hits_file',gene_detector, {})
        refine_contig_neighbours_object.genes = {  'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig2',\
                                            'gene4':'Contig2', 'gene5':'Contig2'}
        refine_contig_neighbours_object.refine_contig_neighbours()
        self.assertEqual(sorted(refine_contig_neighbours_object.refined_neighbouring_contigs), [])
    
    def test_contig_joins_in_middle(self):
        """Tests when an intersection is found in the case where one end of the contig is closest to the middle of a second"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'),\
                                        ('gene4', 'gene5'), ('gene5', 'gene6') ,('gene3', 'gene7'),\
                                        ('gene7', 'gene8')])
                                        
        neighbouring_contigs = [[('Contig1', 'Contig2'),1,['gene3', 'gene7']]]
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/empty_file' )
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/empty_file',gene_detector, {})
        refine_contig_neighbours_object.genes = {  'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1', 'gene4':'Contig1',\
                                            'gene5':'Contig1', 'gene6':'Contig1', 'gene7':'Contig2', 'gene8':'Contig2'}
    
        self.assertEqual(sorted(refine_contig_neighbours_object.refine_contig_neighbours()), [])
    
    def test_short_contig_in_middle(self):
        """Tests the case when there is a contig with only two genes between longer ones, we should find that the connection that needs more iterations is removed"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'),\
                                        ('gene4', 'gene5'), ('gene5', 'gene6'), ('gene6', 'gene7'),\
                                        ('gene7', 'gene8')])
                                        
        neighbouring_contigs = [[('Contig1', 'Contig2'),1,['gene3','gene4']],\
                                [('Contig2', 'Contig3'),1,['gene5', 'gene6']],\
                                [('Contig1', 'Contig3'),2,['gene4', 'gene5']]]
                                
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file' )
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file', gene_detector, {})
        refine_contig_neighbours_object.contigs = { 'Contig1': {'gene1':None, 'gene2':None, 'gene3':None},\
                                                    'Contig2':{'gene4':None, 'gene5':None},\
                                                    'Contig3':{'gene6':None, 'gene7':None, 'gene8':None}}
                                            
        refine_contig_neighbours_object.genes = {   'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1',\
                                                    'gene4':'Contig2', 'gene5':'Contig2', 'gene6':'Contig3',\
                                                    'gene7':'Contig3', 'gene8':'Contig3'}
    
        self.assertEqual(sorted(refine_contig_neighbours_object.refine_contig_neighbours()), sorted([[('Contig1', 'Contig2'),1,['gene3', 'gene4']], [('Contig2', 'Contig3'),1,['gene5', 'gene6']]]))
    
    def test_keep_three_cycle_two_equal_weights(self):
        """Tests that when there are three cycles with all connections having equal weight, these are preserved under refinement"""
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([ ('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'),\
                                        ('gene4', 'gene5'), ('gene6', 'gene7'), ('gene7', 'gene8'),\
                                        ('gene8', 'gene1')])
        
        neighbouring_contigs = [[('Contig1', 'Contig2'), 1, ['gene3', 'gene4']],\
                                [('Contig2', 'Contig3'), 1, ['gene5', 'gene6']],\
                                [('Contig1', 'Contig3'), 1, ['gene8', 'gene1']]]
        
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file')
        
        
        refine_contig_neighbours_object = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file', gene_detector, {})
        
        refine_contig_neighbours_object.contigs = { 'Contig1':{'gene1':None, 'gene2':None, 'gene3':None},\
                                                    'Contig2':{'gene4':None, 'gene5':None}              ,\
                                                    'Contig3':{'gene6':None, 'gene7':None, 'gene8':None}}
                                                    
        refine_contig_neighbours_object.genes = {   'gene1':'Contig1', 'gene2':'Contig1', 'gene3':'Contig1',\
                                                    'gene4':'Contig2', 'gene5':'Contig2', 'gene6':'Contig3',\
                                                    'gene7':'Contig3', 'gene8':'Contig3'}
        
        self.assertEqual(sorted(refine_contig_neighbours_object.refine_contig_neighbours()), sorted(neighbouring_contigs))
    
    
    
    def test_contig_appearances(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1', 'gene2'), ('gene2', 'gene3'), ('gene3', 'gene4'), ('gene4', 'gene5')])
        neighbouring_contigs = [[('Contig1', 'Contig2'),1,['gene3','gene4']]]
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file')
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file', gene_detector, {})
        self.assertDictEqual(refine_contig_neighbours.find_contig_appearances(neighbouring_contigs[0]),{'Contig1':[3, {0:['gene3'], 1:['gene2'], 2:['gene1']}],\
                                                                                                        'Contig2':[2, {0:['gene4'], 1:['gene5']}]})
    
    def test_contig_appearances_multiple_genes_on_one_iteration(self):
        filtered_graph = nx.Graph()
        filtered_graph.add_edges_from([('gene1','gene2'), ('gene2', 'gene3'), ('gene2', 'geneA'), ('geneA', 'gene4'), ('gene4', 'gene5'), ('gene5', 'gene6')])
        neighbouring_contigs =[[('Contig1', 'Contig2'),1,['geneA']], [('Contig2','Contig3'),1,['gene5', 'gene6']]]
        gene_detector = GeneDetector('madansi/tests/data/assembly_4_sequences.fa', 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file')
        refine_contig_neighbours = RefineContigNeighbours(neighbouring_contigs, filtered_graph, 'madansi/tests/data/refine_contig_neighbours_9_blast_hits_file', gene_detector, {})
        self.assertDictEqual(refine_contig_neighbours.find_contig_appearances(neighbouring_contigs[0]),{'Contig1':[3, {1:['gene2'], 2:sorted(['gene1', 'gene3'])}],\
                                                                                                        'Contig2':[2, {1:['gene4'], 2:['gene5']}],\
                                                                                                        'Contig3':[1, {3:['gene6']}]})
        
        