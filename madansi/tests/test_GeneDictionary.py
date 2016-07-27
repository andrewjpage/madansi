import unittest
from madansi.GeneDictionary import GeneDictionary

class TestGeneDictionary(unittest.TestCase):

    def test__init__(self):
        gd = GeneDictionary("ab", "cd")
        self.assertTrue(gd.graphfile)
        self.assertTrue(gd.filteredfile)

    def test_construct_dictionary(self):
        """Tests whether the expected dictionary is the same as the dictionary created and whether the dictionary correctly 
        identifies the presence of some genes and the lack of presence of others"""
        gd = GeneDictionary('madansi/tests/data/graph_4_nodes.dot','madansi/tests/data/gene_present_unittest')
        gene_dict = gd.construct_dictionary()
        expected_dict = {'gene2':(True,'7.23.B265.9.cap3_contig',False), 'gene3':(True, '7.23.B265.9.cap3_contig', True), 'gene1':(True ,'Sample3', False), 'gene4':(False, 'Sample4', True)}
        self.maxdiff = 80
        self.assertDictEqual(gene_dict, expected_dict)
    
    def test_index_file(self):
        gd = GeneDictionary('madansi/tests/data/graph_4_nodes.dot', 'madansi/tests/data/gene_present_unittest')
        index_lines = gd.index_filtered_file()
        expected_dict = {'gene3':1, 'gene2':0, 'gene1':2, 'gene4':3}
        self.assertDictEqual(index_lines, expected_dict)