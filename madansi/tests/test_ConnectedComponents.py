from madansi.ConnectedComponents import ConnectedComponents
from madansi.Component import Component
import networkx as nx
import unittest

class TestConnectedComponents(unittest.TestCase):
    
    
    def test_intitialisation(self):
        components = ConnectedComponents('')
        self.assertEqual(components.components_dictionary, {})
    
    def test_empty_graph(self):
        empty_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/empty_graph.dot'))
        components = ConnectedComponents(empty_graph)
        self.assertEqual(components.components_dictionary, {})
    
    def test_more_complicated_graph(self): #May want to change the asserts equal to look at objects
        test_graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/connected_component_graph.dot'))
        components = ConnectedComponents(test_graph)
        contig_component_1 = Component(0, sorted(['contig4', 'contig5', 'contig6', 'contig7']), sorted(['contig4', 'contig7']),\
                                         sorted([('contig4', 'contig5'), ('contig5', 'contig6'), ('contig6', 'contig7')]))
        contig_component_2 = Component(1, sorted(['contig2', 'contig3']), sorted(['contig2', 'contig3']), [('contig2', 'contig3')])
        contig_component_3 = Component(2, ['contig1'], ['contig1'], [])
        component_dictionary = components.connected_components()
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