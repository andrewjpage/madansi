from madansi.ConnectedComponents import ConnectedComponents
from madansi.Component import Component
import networkx as nx
import unittest

class TestConnectedComponents(unittest.TestCase):
    
    def test_empty_graph(self):
        empty_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/empty_graph.dot'))
        components = ConnectedComponents(empty_graph , empty_graph)
        self.assertEqual(components.create_component_dictionary(), {})
    
    def test_component_dictionary(self): #May want to change the asserts equal to look at objects
        test_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/connected_component_graph.dot'))
        components = ConnectedComponents(test_graph , '')
        contig_component_1 = Component(0, sorted(['contig4', 'contig5', 'contig6', 'contig7']), sorted(['contig4', 'contig7']),\
                                         sorted([('contig4', 'contig5'), ('contig5', 'contig6'), ('contig6', 'contig7')]))
        contig_component_2 = Component(1, sorted(['contig2', 'contig3']), sorted(['contig2', 'contig3']), [('contig2', 'contig3')])
        contig_component_3 = Component(2, ['contig1'], ['contig1'], [])
        component_dictionary = components.create_component_dictionary()
        self.assertEqual(sorted(component_dictionary.keys()), sorted([0,1,2]))
        self.assertEqual(sorted(component_dictionary[0].contigs), sorted(['contig4', 'contig5', 'contig6', 'contig7']) )
        self.assertEqual(sorted(component_dictionary[0].ends), sorted(['contig4', 'contig7']))
        self.assertEqual(sorted(component_dictionary[0].edges), sorted([tuple(sorted(('contig4', 'contig5'))), tuple(sorted(('contig5', 'contig6'))), tuple(sorted(('contig6', 'contig7')))]))
        self.assertEqual(sorted(component_dictionary[1].contigs), sorted(['contig2', 'contig3']))
        self.assertEqual(sorted(component_dictionary[1].ends), sorted(['contig2', 'contig3']))
        self.assertEqual(sorted(component_dictionary[1].edges), [tuple(sorted(('contig2', 'contig3')))])
        self.assertEqual(sorted(component_dictionary[2].contigs), ['contig1'])
        self.assertEqual(sorted(component_dictionary[2].ends), ['contig1'])
        self.assertEqual(sorted(component_dictionary[2].edges), [])
    
    def test_contigs_to_components(self):
        test_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/connected_component_graph.dot'))
        components = ConnectedComponents(test_graph, test_graph)
        expected_dictionary = {'contig1':2, 'contig2':1,'contig3':1, 'contig4':0, 'contig5':0, 'contig6':0, 'contig7':0}
        self.assertDictEqual(components.create_contigs_to_components_dictionary(), expected_dictionary)
    
    def test_contigs_to_components_additional_contigs(self):
        refined_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/connected_component_graph.dot'))
        unrefined_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/connected_component_graph_extra_contigs.dot'))
        components = ConnectedComponents(refined_graph, unrefined_graph)
        expected_dictionary = {'contig1':2, 'contig2':1,'contig3':1, 'contig4':0, 'contig5':0, 'contig6':0, 'contig7':0, 'contig8':None, 'contig9':None}
        self.assertDictEqual(components.create_contigs_to_components_dictionary(), expected_dictionary)
        