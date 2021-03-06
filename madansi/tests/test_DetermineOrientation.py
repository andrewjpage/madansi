import unittest
from madansi.DetermineOrientation import DetermineOrientation
import networkx as nx

class TestDetermineOrientation(unittest.TestCase):
    
    def test_determine_orientation(self):
        """Tests the basic case with tuples and start, end and middle contigs"""
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3'), ('Contig3', 'Contig4')])
        contig_ends = {'Contig1':{'Contig2': (1,100)}, 'Contig2':{'Contig1':(29, 102), 'Contig3':(307, 240)}, 'Contig3':{'Contig2':(1000,891), 'Contig4':(201, 401)}, 'Contig4':{'Contig3':(300, 49)}}
        sequences = {'Contig1':['','', 100], 'Contig2':['','', 400], 'Contig3':['','',1200], 'Contig4':['','',400]}
        
        self.assertEqual(DetermineOrientation().determine_orientation('Contig3', contig_ends, ['Contig1', 'Contig2'], {}, graph, False, sequences), -1)
        self.assertEqual(DetermineOrientation().determine_orientation_end_contig('Contig4', contig_ends,['Contig1', 'Contig2', 'Contig3'] , {}, graph, False, sequences), -1)
        self.assertEqual(DetermineOrientation().determine_orientation('Contig4', contig_ends,['Contig1', 'Contig2', 'Contig3'] , {}, graph, False, sequences), -1)
        self.assertEqual(DetermineOrientation().determine_orientation('Contig2', contig_ends, ['Contig1'], {}, graph, False, sequences), 1)
        self.assertEqual(DetermineOrientation().determine_orientation_start_contig('Contig1', contig_ends, [], {}, graph, False, sequences), -1)
        
    
    def test_determine_orientation_None_second_entry(self):
        """Tests when only one gene is found when walking out"""
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3')])
        contig_ends = {'Contig1':{'Contig2':(1,None)}, 'Contig2':{'Contig1':(29,102), 'Contig3':(307,None)}, 'Contig3':{'Contig2':(400, None)}}
        sequences = {'Contig1':['', '',500], 'Contig2':['','',971], 'Contig3':['', '', 500]}
        
        self.assertEqual(DetermineOrientation().determine_orientation_start_contig('Contig1', contig_ends, [], {}, graph, False, sequences), -1)
        self.assertEqual(DetermineOrientation().determine_orientation('Contig2', contig_ends, ['Contig1'], {}, graph, False, sequences), 1)
        self.assertEqual(DetermineOrientation().determine_orientation('Contig3', contig_ends, [], {}, graph, False, sequences), -1)
        self.assertEqual(DetermineOrientation().determine_orientation_end_contig('Contig3', contig_ends, [], {}, graph, False, sequences), -1)
    
    def test_contig_orientation_f_f(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(490,3)}, 'Contig2': {'Contig1':(1,300)}}
        sequences = {'Contig1':['', '', 500], 'Contig2':['','',300]}
        
        self.assertEqual(DetermineOrientation().determine_orientation_start_contig('Contig1', contig_ends, [], {}, graph, False, sequences), 1)
        self.assertEqual(DetermineOrientation().determine_orientation_end_contig('Contig2', contig_ends, ['Contig1'], {}, graph, False, sequences), 1)
    
    def test_contig_orientation_f_r(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(490,3)}, 'Contig2': {'Contig1':(300,1)}}
        sequences = {'Contig1':['', '', 500], 'Contig2':['','',300]}
        
        self.assertEqual(DetermineOrientation().determine_orientation_start_contig('Contig1', contig_ends, [], {}, graph, False, sequences), 1)
        self.assertEqual(DetermineOrientation().determine_orientation_end_contig('Contig2', contig_ends, ['Contig1'], {}, graph, False, sequences), -1)
    
    def test_contig_orientation_r_r(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(3,490)}, 'Contig2': {'Contig1':(300,1)}}
        sequences = {'Contig1':['', '', 500], 'Contig2':['','',300]}
       
        self.assertEqual(DetermineOrientation().determine_orientation_start_contig('Contig1', contig_ends, [], {}, graph, False, sequences), -1)
        self.assertEqual(DetermineOrientation().determine_orientation_end_contig('Contig2', contig_ends, ['Contig1'], {}, graph, False, sequences), -1)
    
    def test_contig_orientation_r_f(self):
        graph = nx.Graph()
        graph.add_edge('Contig1', 'Contig2')
        contig_ends = {'Contig1':{'Contig2':(3,490)}, 'Contig2': {'Contig1':(1,300)}}
        sequences = {'Contig1':['', '', 500], 'Contig2':['','',300]}
        
        self.assertEqual(DetermineOrientation().determine_orientation_start_contig('Contig1', contig_ends, [], {}, graph, False, sequences), -1)
        self.assertEqual(DetermineOrientation().determine_orientation_end_contig('Contig2', contig_ends, ['Contig1'], {}, graph, False, sequences), 1)
    
    def test_sections_contig_reordered(self):
        """Tests if the contig splits so that one of the original ends of the contig now lies in the middle"""
        graph = nx.Graph()
        graph.add_edges_from([('Contig1', 'Contig2'), ('Contig2', 'Contig3')])
        contig_ends = {'Contig1': {'Contig2':(1000,200)}, 'Contig2':{'Contig1':(3564,2937), 'Contig3':(9099,9491)}, 'Contig3':{'Contig2':(10, 918)}}
        
        self.assertEqual(DetermineOrientation().determine_orientation_degree_2('Contig2', contig_ends, ['Contig3'], {}, graph, False, {'Contig1':['','',1000], 'Contig2':['','', 40000], 'Contig3':['','','',1000]}), 1)
        
        
        
