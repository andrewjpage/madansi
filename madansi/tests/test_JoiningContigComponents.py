import unittest
import networkx as nx
from madansi.JoiningContigComponents import JoiningContigComponents

class TestJoiningContigComponents(unittest.TestCase):
    
    def test_list_ends_of_components(self):
        """Tests that all of the ends of the contig components are found"""
        refined_graph = nx.Graph()
        refined_graph.add_edges_from([('contig1', 'contig2'), ('contig5', 'contig6')])
        unrefined_graph = nx.Graph()
        unrefined_graph.add_edges_from([('contig1', 'contig2'), ('contig2', 'contig3'),\
                                        ('contig3', 'contig4'), ('contig4', 'contig5'),\
                                        ('contig5', 'contig6')])
        expected_ends_of_components = ['contig1', 'contig2', 'contig5', 'contig6']
        
        join_contig_components = JoiningContigComponents(refined_graph, unrefined_graph)
        ends_of_components = join_contig_components.list_ends_of_components()
        self.assertEqual(sorted(expected_ends_of_components), sorted(ends_of_components))
        
    def test_direct_joining_of_components(self):
        """Tests the case when two contig components are adjacent"""
        refined_graph = nx.Graph()
        refined_graph.add_edges_from([('contig1', 'contig2'), ('contig3', 'contig4')])
        unrefined_graph = nx.Graph()
        unrefined_graph.add_edges_from([('contig1', 'contig2'), ('contig2', 'contig3'), ('contig3', 'contig4')])
        
        join_contig_components  = JoiningContigComponents(refined_graph, unrefined_graph)
        ordered_contig_graph = join_contig_components.add_edges()
        self.assertTrue(nx.is_isomorphic(ordered_contig_graph, unrefined_graph))
    
    def test_add_edges_from_ends_with_degree_two(self):
        """Tests the case when there are single contigs adjacent to a component"""
        refined_graph = nx.Graph()
        refined_graph.add_edges_from([('contig1', 'contig2'), ('contig5', 'contig6')])
        unrefined_graph = nx.Graph()
        unrefined_graph.add_edges_from([('contig1', 'contig2'), ('contig2', 'contig3'),\
                                        ('contig3', 'contig4'), ('contig4', 'contig5'),\
                                        ('contig5', 'contig6')])
        expected_ordered_contig_graph = nx.Graph()
        expected_ordered_contig_graph.add_edges_from([  ('contig1', 'contig2'), ('contig2', 'contig3'),\
                                                        ('contig4', 'contig5'), ('contig5', 'contig6')])
        
        join_contig_components = JoiningContigComponents(refined_graph, unrefined_graph)
        ordered_contig_graph = join_contig_components.add_edges()
        self.assertTrue(nx.is_isomorphic(ordered_contig_graph, expected_ordered_contig_graph))
    
    def test_end_of_component_meets_middle_of_component(self):
        """Tests when one contig component meets the middle of a second and this should not make a change to the graph"""
        refined_graph = nx.Graph()
        refined_graph.add_edges_from([('contig1', 'contig2'), ('contig2', 'contig3'), ('contig4', 'contig5')])
        unrefined_graph = nx.Graph()
        unrefined_graph.add_edges_from([('contig1', 'contig2'), ('contig2', 'contig3'),\
                                        ('contig2', 'contig4'), ('contig4', 'contig5')])
        
        join_contig_components  = JoiningContigComponents(refined_graph, unrefined_graph)
        ordered_contig_graph = join_contig_components.add_edges()
        self.assertTrue(nx.is_isomorphic(ordered_contig_graph, refined_graph))