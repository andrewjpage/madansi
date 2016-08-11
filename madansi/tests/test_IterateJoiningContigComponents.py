import unittest
import networkx as nx
from madansi.IterateJoiningContigComponents import IterateJoiningContigComponents

class TestIterateJoiningContigComponents(unittest.TestCase):
    
    def test_empty_graph(self):
        """Tests when both the unrefined and refined graphs are empty"""
        test_case = IterateJoiningContigComponents(nx.Graph(),'')
        ordered_contig_graph = test_case.iterate_joining(nx.Graph())
        self.assertTrue(nx.is_isomorphic(nx.Graph(), ordered_contig_graph))
    
    def test_empty_initial_refined_graph(self):
        """Tests when the refined graph is initially empty but the unrefined graph is not"""
        unrefined_graph = nx.Graph()
        unrefined_graph.add_edge('contig1', 'contig2')
        test_case = IterateJoiningContigComponents(unrefined_graph, '')
        ordered_contig_graph = test_case.iterate_joining(nx.Graph())
        self.assertTrue(nx.is_isomorphic(nx.Graph(), ordered_contig_graph))
    
    def test_one_iteration(self):
        """Tests when only one iteration is required"""
        unrefined_graph = nx.Graph()
        unrefined_graph.add_edges_from([('contig1', 'contig2') , ('contig2', 'contig3')])
        initial_refined_graph = nx.Graph()
        initial_refined_graph.add_edge('contig1', 'contig2')
        
        test_case = IterateJoiningContigComponents(unrefined_graph, '')
        ordered_contig_graph = test_case.iterate_joining(initial_refined_graph)
        self.assertTrue(nx.is_isomorphic(ordered_contig_graph, unrefined_graph))
    
    def test_three_iterations(self):
        """Tests when three iterations are needed until the process terminates"""
        unrefined_graph = nx.Graph()
        unrefined_graph.add_edges_from([('contig1', 'contig2'), ('contig2', 'contig3'),\
                                        ('contig3', 'contig4'), ('contig4', 'contig5')])
        initial_refined_graph = nx.Graph()
        initial_refined_graph.add_edge('contig1', 'contig2')
        
        test_case = IterateJoiningContigComponents(unrefined_graph, '')
        ordered_contig_graph = test_case.iterate_joining(initial_refined_graph)
        self.assertTrue(nx.is_isomorphic(ordered_contig_graph, unrefined_graph))
    
    def test_finishes_different_to_unrefined_graph(self):
        unrefined_graph = nx.Graph()
        initial_refined_graph = nx.Graph()
        final_refined_graph = nx.Graph()
        
        unrefined_graph.add_edges_from([('contig1', 'contig2'), ('contig2', 'contig3'),\
                                        ('contig3', 'contig4'), ('contig4', 'contig5'),\
                                        ('contig4', 'contig6')])
        
        initial_refined_graph.add_edge('contig1', 'contig2')
        initial_refined_graph.add_node('contig6')
        
        final_refined_graph.add_edges_from([('contig1','contig2'), ('contig2', 'contig3'),\
                                            ('contig3','contig4'), ('contig4', 'contig6')])
        
        test_case = IterateJoiningContigComponents(unrefined_graph, '')
        ordered_contig_graph = test_case.iterate_joining(initial_refined_graph)
        print(ordered_contig_graph.edges())
        self.assertTrue(nx.is_isomorphic(ordered_contig_graph, final_refined_graph))
        
        