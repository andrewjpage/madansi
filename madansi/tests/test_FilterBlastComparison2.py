import unittest
from madansi.FilterBlastComparison2 import FilterBlastComparison2
import os
import filecmp

class TestFilterBlastComparison(unittest.TestCase):
    
    def test__init__FilterBlastComaparison(self):
        fbc = FilterBlastComparison2('ab','cd', 20.0, 5, bitscore=23.4)
        self.assertTrue(fbc.comparisonfile)
        self.assertTrue(fbc.filteredfile)
        self.assertIsInstance(fbc.percentidentity, float)
        self.assertIsInstance(fbc.alignmentlength, int)
        self.assertIsInstance(fbc.mismatches ,int)
        self.assertIsInstance(fbc.gapopenings, int)
        self.assertIsInstance(fbc.evalue, float)
        self.assertIsInstance(fbc.bitscore, float)
        
    def test_filter(self):
        output_filter_percent_identity = 'output'
        fbc = FilterBlastComparison2('madansi/tests/data/blast_unittest.m8', output_filter_percent_identity, percentidentity=99)
        fbc.filter_percent_identity()
        self.assertTrue(filecmp.cmp(output_filter_percent_identity, 'madansi/tests/data/expected_filter_percent_identity'))
        self.clean_up()
        
    def test_filter_alignment_length(self):
        output_filter_alignment_length = 'output'
        fbc = FilterBlastComparison2('madansi/tests/data/blast_unittest.m8', output_filter_alignment_length, alignmentlength=230)
        fbc.filter_alignment_length()
        self.assertTrue(filecmp.cmp(output_filter_alignment_length, 'madansi/tests/data/expected_filter_alignment_length'))
        self.clean_up()
        
    def test_filter_mismatches(self):
        output_filter_mismatches = 'output'
        fbc = FilterBlastComparison2('madansi/tests/data/blast_unittest.m8', output_filter_mismatches, mismatches=2.8)
        fbc.filter_mismatches()
        self.assertTrue(filecmp.cmp(output_filter_mismatches, 'madansi/tests/data/expected_filter_mismatches'))
        self.clean_up()
        
        
#    def test_filter_gapopenings(self):
#        output_filter_gapopenings = 'output'
#        fbc = FilterBlastComparison2('madansi/tests/data/blast_unittest.m8', output_filter_gapopenings, gapopenings=2)
#        fbc.filter_gapopenings()
#        self.assertTrue(filecmp.cmp(output_filter_gapopenings, 'madansi/tests/data/expected_filter_gapopenings'))
        
    def test_filter_evalue(self):
        output_filter_evalue= 'output'
        fbc = FilterBlastComparison2('madansi/tests/data/blast_unittest.m8', output_filter_evalue, evalue=1e-100)
        fbc.filter_evalue()
        self.assertTrue(filecmp.cmp(output_filter_evalue, 'madansi/tests/data/expected_filter_evalue'))
        self.clean_up()
        
    def test_filter_bit_score(self):
        output_filter_bit_score = 'output'
        fbc = FilterBlastComparison2('madansi/tests/data/blast_unittest.m8', output_filter_bit_score, bitscore=400)
        fbc.filter_bit_score()
        self.assertTrue(filecmp.cmp(output_filter_bit_score, 'madansi/tests/data/expected_filter_bit_score'))
        self.clean_up()

    def clean_up(self):
        if os.path.exists('output'):
            os.remove('output')
    
    
    