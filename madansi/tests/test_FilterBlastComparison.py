import unittest
from madansi.FilterBlastComparison import FilterBlastComparison
import os
import filecmp

class TestFilterBlastComparison(unittest.TestCase):
    
    def test__init__FilterBlastComaparison(self):
        fbc = FilterBlastComparison('ab','cd', percentidentity=20.0, alignmentlength=5, gapopenings=0, mismatches=2, evalue=17.3, bitscore=23.4)
        self.assertEqual(fbc.comparisonfile,'ab')
        self.assertEqual(fbc.filteredfile, 'cd')
        self.assertEqual(fbc.percentidentity, 20.0)
        self.assertEqual(fbc.alignmentlength, 5)
        self.assertEqual(fbc.mismatches ,2)
        self.assertEqual(fbc.gapopenings, 0)
        self.assertEqual(fbc.evalue, 17.3)
        self.assertEqual(fbc.bitscore, 23.4)
        
    def test_filter(self):
        """Testing using multiple filters on the data and values given without name"""
        output_filter = 'output'
        fbc = FilterBlastComparison('madansi/tests/data/blast_unittest.m8', output_filter, alignmentlength=230, bitscore=350)
        fbc.filter()
        self.assertTrue(filecmp.cmp(output_filter, 'madansi/tests/data/expected_filter'))
        os.unlink(output_filter)
        
        output_filter = 'output'
        fbc = FilterBlastComparison('madansi/tests/data/blast_unittest.m8', output_filter, 99)
        fbc.filter()
        self.assertTrue(filecmp.cmp(output_filter, 'madansi/tests/data/expected_filter_test_2'))
        os.unlink(output_filter)
    
    def test_no_filters(self):
        """No optional inputs given outputs the original data"""
        output_filter = 'output'
        fbc = FilterBlastComparison('madansi/tests/data/blast_unittest.m8', output_filter)
        fbc.filter()
        self.assertTrue(filecmp.cmp(output_filter, 'madansi/tests/data/blast_unittest.m8'))
        os.unlink(output_filter)
        
    def tests_invalid_percentage_identity(self):
        """Tests a percentage identity outside the valid range"""
        fbc = FilterBlastComparison('madansi/tests/data/blast_unittest.m8', 'output', percentidentity=103.0)
        self.assertRaises(ValueError)
