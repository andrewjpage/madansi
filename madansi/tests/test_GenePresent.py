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
        expected_dict = {'Contig2':False, 'Contig3':True, '7.23.B265.14.cap3_contig':True , '7.23.B265.9.cap3_contig':False}
        self.assertDictEqual(gene_present_dict, expected_dict)

        self.assertEqual(True, gene_present_dict['7.23.B265.14.cap3_contig'])
        self.assertEqual(gene_present_dict['Contig2'], False)

    def test_list_genes_present(self):
        """Tests the list of genes present generated from the dictionary"""
        gene_present_dict = {'gene1':True, 'gene2':True, 'gene3':False, 'gene4':True, 'gene5':False}
        expected_list = ['gene1', 'gene2', 'gene4']
        new_list = list_genes_present(gene_present_dict)
        self.assertCountEqual(expected_list, new_list)