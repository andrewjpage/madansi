import unittest
from madansi.Component import Component

class TestComponent(unittest.TestCase):
    def test_initialisation(self):
        component = Component(0, ['contig1', 'contig2', 'contig3', 'contig4'], ['contig1', 'contig4'], [('contig1', 'contig2'), ('contig2', 'contig3'), ('contig3', 'contig4')])
        self.assertEqual(component.component_index, 0)
        self.assertEqual(sorted(component.contigs), sorted(['contig1', 'contig2', 'contig3', 'contig4']))
        self.assertEqual(sorted(component.ends), sorted(['contig1', 'contig4']))
        self.assertEqual(sorted(component.edges), sorted([('contig1', 'contig2'), ('contig2', 'contig3'), ('contig3', 'contig4')]))