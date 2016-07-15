import unittest
from madansi.FilterBlastComparison import FilterBlastComparison
import os
import filecmp

class TestFilterBlastComparison(unittest.TestCase):
    
    def test__init__FilterBlastComaparison(self):
        fbc = FilterBlastComparison('ab','cd', 20.0, 5, bitscore=23.4)
        self.assertTrue(fbc.comparisonfile)
        self.assertTrue(fbc.filteredfile)
        self.assertIsInstance(fbc.percentidentity, float)
        self.assertIsInstance(fbc.alignmentlength, int)
        self.assertIsInstance(fbc.mismatches ,int)
        self.assertIsInstance(fbc.gapopenings, int)
        self.assertIsInstance(fbc.evalue, float)
        self.assertIsInstance(fbc.bitscore, float)
        
    def test_filter(self):
        output_filter = 'output'
        fbc = FilterBlastComparison('madansi/tests/data/blast_unittest.m8', output_filter, alignmentlength=230, bitscore=350)
        fbc.filter()
        self.assertTrue(filecmp.cmp(output_filter, 'madansi/tests/data/expected_filter'))
        self.clean_up()
            
    def clean_up(self):
        if os.path.exists('output'):
            os.remove('output')