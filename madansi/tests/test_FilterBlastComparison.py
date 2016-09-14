import unittest
from madansi.FilterBlastComparison import FilterBlastComparison
import os
import filecmp
import shutil
import tempfile

class TestFilterBlastComparison(unittest.TestCase):
    
    def test__init__FilterBlastComaparison(self):
        """Tests initialisation of all parameters, including both those specified and not"""
        temp_dir = tempfile.mkdtemp()
        fbc = FilterBlastComparison('ab',temp_dir, percent_identity=20.0, alignment_length=5, gap_openings=0, mismatches=2, evalue=17.3, bit_score=23.4)
        self.assertEqual(fbc.input_blast_file,'ab')
        self.assertTrue(fbc.filtered_blast_output)
        self.assertEqual(fbc.percent_identity, 20.0)
        self.assertEqual(fbc.alignment_length, 5)
        self.assertEqual(fbc.mismatches ,2)
        self.assertEqual(fbc.gap_openings, 0)
        self.assertEqual(fbc.evalue, 17.3)
        self.assertEqual(fbc.bit_score, 23.4)       
        shutil.rmtree(temp_dir)
    
    def test_filter_duplicate_genes(self):
        """Tests that only the genes that appear on one contig are outputted with no explicit filters on the other values"""
        temp_dir = tempfile.mkdtemp()
        fbc = FilterBlastComparison('madansi/tests/data/blast_unittest.m8', temp_dir)
        fbc.filter()
        self.assertTrue(filecmp.cmp(fbc.filtered_blast_output.name, 'madansi/tests/data/expected_filter_remove_duplicates', shallow=False))
        shutil.rmtree(temp_dir)
    
    def test_filter_no_duplicate_genes(self):
        """Tests the case when there are no duplicate genes and no explicit filters given so will return the original list"""
        temp_dir = tempfile.mkdtemp()
        fbc = FilterBlastComparison('madansi/tests/data/test_blast_hits', temp_dir)
        fbc.filter()
        self.assertTrue(filecmp.cmp(fbc.filtered_blast_output.name, 'madansi/tests/data/test_blast_hits', shallow=False))
        shutil.rmtree(temp_dir)
    
    def test_filter_by_values(self):
        """Testing using multiple filters on the data"""
        temp_dir = tempfile.mkdtemp()
        fbc = FilterBlastComparison('madansi/tests/data/test_blast_hits', temp_dir, 99, alignment_length = 230, bit_score = 350)
        fbc.filter()
        self.assertTrue(filecmp.cmp(fbc.filtered_blast_output.name, 'madansi/tests/data/expected_filter_values', shallow=False))
        shutil.rmtree(temp_dir)
        
    def test_duplicate_genes_and_filters(self):
        """Tests that the correct filtered file is obtained when both types of filter have to be used"""
        temp_dir = tempfile.mkdtemp()
        fbc = FilterBlastComparison('madansi/tests/data/blast_unittest.m8', temp_dir, bit_score = 500)
        fbc.filter()
        self.assertTrue(filecmp.cmp(fbc.filtered_blast_output.name, 'madansi/tests/data/empty_file', shallow=False))
        shutil.rmtree(temp_dir)
        
    def tests_invalid_percentage_identity(self):
        """Tests a percentage identity outside the valid range"""
        temp_dir = tempfile.mkdtemp()
        FilterBlastComparison('madansi/tests/data/blast_unittest.m8', temp_dir, percent_identity=103.0)
        self.assertRaises(ValueError)
        shutil.rmtree(temp_dir)
