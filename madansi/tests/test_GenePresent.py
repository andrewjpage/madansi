import unittest
from madansi.GenePresent import GenePresent, list_genes_present

class TestGenePresent(unittest.TestCase):

    def test__init__(self):
        gp = GenePresent("ab")
        self.assertTrue(gp.filteredfile)

    def test_construct_dictionary(self):
        """Tests whether the expected dictionary is the same as the dictionary created and whether the dictionary correctly 
        identifies the presence of some genes and the lack of presence of others"""
        gp = GenePresent('madansi/tests/data/gene_present_unittest')
        gene_present_dict = gp.construct_dictionary()
        expected_dict = {'gene2':True, 'gene3':True, 'gene1':True , 'gene4':False}
        self.assertDictEqual(gene_present_dict, expected_dict)

        self.assertEqual(True, gene_present_dict['gene1'])
        self.assertEqual(gene_present_dict['gene4'], False)

    def test_list_genes_present(self):
        """Tests the list of genes present generated from the dictionary"""
        gene_present_dict = {'gene1':True, 'gene2':True, 'gene3':False, 'gene4':True, 'gene5':False}
        expected_list = ['gene1', 'gene2', 'gene4']
        new_list = list_genes_present(gene_present_dict)
        self.assertCountEqual(expected_list, new_list)
    
    def test_index_file(self):
        gp = GenePresent('madansi/tests/data/gene_present_unittest')
        index_lines = gp.index_filtered_file()
        expected_dict = {'gene3':1, 'gene2':0, 'gene1':2, 'gene4':3}
        self.assertDictEqual(index_lines, expected_dict)