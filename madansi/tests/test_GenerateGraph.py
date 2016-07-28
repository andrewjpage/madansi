import unittest
from madansi.GenerateGraph import GenerateGraph
import networkx as nx
import filecmp
import os

class TestGenerateGraph(unittest.TestCase):
    
    def test__init__(self):
        gg = GenerateGraph('ab','cd','ef','gh')
        self.assertTrue(gg.graphfile)
        self.assertTrue(gg.filteredfile)
        self.assertTrue(gg.outputgraphfile)
        self.assertTrue(gg.unused_sequence_file)
    
    def test_open_graph_file(self):
        gg = GenerateGraph('madansi/tests/data/graph_3_nodes.dot', 'cd', 'ef', 'gh')
        g = gg.open_graph_file()
        self.assertTrue(g)
    
    def test_create_dictionaries(self):
        gg = GenerateGraph('ab', 'madansi/tests/data/gene_present_unittest', 'ef', 'gh')
        index_file = gg.create_index_file()
        expected_dict = {'gene2':0, 'gene3':1, 'gene1':2, 'gene4':3}
        self.assertDictEqual(expected_dict,index_file)
        
        gene_present = gg.create_gene_present_dict()
        expected_dict = {'gene2': True, 'gene3':True, 'gene1':True, 'gene4':False}
        self.assertDictEqual(gene_present, expected_dict)
        
        gene_sequence = gg.create_gene_sequence_dict()
        expected_dict = {'gene2':'7.23.B265.9.cap3_contig', 'gene3':'7.23.B265.9.cap3_contig', 'gene1':'Sample3', 'gene4':'Sample4'}
        self.assertDictEqual(gene_sequence, expected_dict)
        
        gene_orientation = gg.create_gene_orientation_dict()
        expected_dict = {'gene2': False, 'gene3': True, 'gene1': False, 'gene4':True}
        self.assertDictEqual(gene_orientation, expected_dict)
    
    def test_ends_of_sequence(self):
        gg = GenerateGraph('madansi/tests/data/graph_5_nodes_2.dot', 'madansi/tests/data/filtered_data_5_contigs_2', 'ef', 'gh')
        gene_present = {'gene2': True, 'gene3':True, 'gene1':True, 'gene4':True, 'gene5': True}
        gene_sequence = {'gene2':'7.23.B265.9.cap3_contig', 'gene3':'7.23.B265.9.cap3_contig', 'gene1':'7.23.B265.9.cap3_contig', 'gene4':'7.23.B265.9.cap3_contig', 'gene5':'7.23.B265.9.cap3_contig'}
        graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_5_nodes_2.dot'))
        
        end_list = gg.ends_of_sequence('gene1', graph , gene_present, gene_sequence )
        expected_list = ['gene1', 'gene5']
        self.assertCountEqual(end_list, expected_list)
        
        end_list = gg.ends_of_sequence('gene3', graph, gene_present, gene_sequence)
        self.assertCountEqual(end_list, expected_list)
        
    def test_closest_gene(self):
        #Linear:   
        gg = GenerateGraph('madansi/tests/data/graph_order_sequences.dot', 'madansi/tests/data/filtered_data_test_order_sequences', 'ef', 'gh')
        gene_present = {'gene1':True, 'gene2':True, 'gene3':True, 'gene4':False, 'gene5':False, 'gene6':True, 'gene7':True, 'gene8':True, 'gene9':False, 'gene10':False}
        gene_sequence = {'gene1':'Sequence1', 'gene2':'Sequence1', 'gene3':'Sequence1', 'gene4':'Sequence1', 'gene5':'Sequence1', 'gene6':'Sequence2', 'gene7':'Sequence2', 'gene8':'Sequence2', 'gene9':'Sequence2', 'gene10':'Sequence2'}
        graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_order_sequences.dot'))
        
        gene = gg.closest_gene('gene1', graph , gene_present, gene_sequence)
        self.assertEqual(gene,'gene6')
        
        gene = gg.closest_gene('gene3', graph , gene_present, gene_sequence)
        self.assertEqual(gene,'gene6')
        
        gene = gg.closest_gene('gene8', graph, gene_present, gene_sequence)
        self.assertEqual(gene,'gene3')
        
        gene = gg.closest_gene('gene6', graph, gene_present, gene_sequence)
        self.assertEqual(gene,'gene3')
        
        #Cycle:
        gg = GenerateGraph('madansi/tests/data/graph_order_cycle.dot', 'madansi/tests/data/filtered_data_test_order_sequences', 'ef', 'gh')
        graph = nx.Graph(nx.drawing.nx_pydot.read_dot('madansi/tests/data/graph_order_cycle.dot'))
        
        gene = gg.closest_gene('gene1', graph, gene_present, gene_sequence)
        self.assertEqual(gene, 'gene8')
        
        gene = gg.closest_gene('gene8', graph, gene_present, gene_sequence)
        self.assertEqual(gene, 'gene1')
        
        gene = gg.closest_gene('gene6', graph, gene_present, gene_sequence)
        self.assertEqual(gene, 'gene3')
        
        gene = gg.closest_gene('gene6', graph, gene_present, gene_sequence)
        self.assertEqual(gene, 'gene3')
    
    
    def test_one_sequence(self):
        """Tests one sequence with no genes- should return the sequence"""
        output_file = 'output.dot'
        unused_sequences = 'outputfile'
        expected_graph = 'madansi/tests/data/GenerateGraph_testdata/expected_one_sequence.dot'
        input_file = 'madansi/tests/data/GenerateGraph_testdata/one_sequence.dot'
        test_data = 'madansi/tests/data/GenerateGraph_testdata/GenerateGraph_test_data'
        expected_unused = 'madansi/tests/data/GenerateGraph_testdata/expected_unused_sequences_one_sequence'
        gg = GenerateGraph(input_file ,test_data , output_file, unused_sequences)
        gg.generate_graph()
        output_list = gg.generate_graph()
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(output_file))
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph))
        self.assertTrue(nx.is_isomorphic(g,h))
        self.assertCountEqual([],output_list)
        os.unlink(output_file)
        os.unlink(unused_sequences)
    
    def test_two_sequences_no_genes(self):
        """Should return an empty graph with both sequences in the unused sequences list"""
        output_file = 'output.dot'
        unused_sequences = 'outputfile'
        expected_graph = 'madansi/tests/data/GenerateGraph_testdata/expected_two_sequences_no_genes.dot'
        input_file= 'madansi/tests/data/GenerateGraph_testdata/two_sequences_no_genes.dot'
        test_data = 'madansi/tests/data/GenerateGraph_testdata/GenerateGraph_test_data'
        gg = GenerateGraph(input_file , test_data, output_file, unused_sequences)
        output_list = gg.generate_graph()
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(output_file))
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph))
        self.assertCountEqual(output_list,['Sequence1', 'Sequence2'])
        self.assertTrue(nx.is_isomorphic(g,h))
        os.unlink(output_file)
        os.unlink(unused_sequences)
    
    def test_two_sequences_one_gene(self):
        """Should return a graph with just the sequence with the whole of the sequence with a gene present and the other sequence should go into the unused sequences list"""
        output_file = 'output.dot'
        unused_sequences = 'outputfile'
        expected_graph = 'madansi/tests/data/GenerateGraph_testdata/expected_two_sequences_one_gene.dot'
        input_file= 'madansi/tests/data/GenerateGraph_testdata/two_sequences_one_gene.dot'
        test_data = 'madansi/tests/data/GenerateGraph_testdata/GenerateGraph_test_data'
        gg = GenerateGraph(input_file , test_data, output_file, unused_sequences)
        output_list = gg.generate_graph()
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(output_file))
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph))
        self.assertTrue(output_list == ['Sequence2'])
        self.assertTrue(nx.is_isomorphic(g,h))
        os.unlink(output_file)
        os.unlink(unused_sequences)
    
   
    def test_two_sequences_two_genes(self):
        """Given two sequences both with genes on will check that the linear combination of the two is given"""
        output_file = 'output.dot'
        unused_sequences = 'outputfile'
        expected_graph = 'madansi/tests/data/GenerateGraph_testdata/expected_two_sequences_two_genes_1.dot'
        input_file= 'madansi/tests/data/GenerateGraph_testdata/two_sequences_two_genes_1.dot'
        test_data = 'madansi/tests/data/GenerateGraph_testdata/GenerateGraph_test_data'
        gg = GenerateGraph(input_file , test_data, output_file, unused_sequences)
        output_list = gg.generate_graph()
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(output_file))
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph))
        self.assertTrue(output_list == [])
        self.assertTrue(nx.is_isomorphic(g,h))
        os.unlink(output_file)
        os.unlink(unused_sequences)
        
    def test_neighbour_nodes_not_in_file(self):
        """Given a graph with present genes surrounded by genes not present in the file, will test first that the correct dictionary of closest pairs is given """
        output_file = 'output.dot'
        unused_sequences = 'outputfile'
        expected_graph_1 = 'madansi/tests/data/GenerateGraph_testdata/expected_neighbor_nodes_not_in_file_1.dot'
        input_file = 'madansi/tests/data/GenerateGraph_testdata/neighbor_nodes_not_in_file.dot'
        test_data = 'madansi/tests/data/GenerateGraph_testdata/GenerateGraph_test_data'        
        gg = GenerateGraph(input_file, test_data, output_file, unused_sequences)
        
        gene_present = {'gene1':True, 'gene2':True, 'gene3':True, 'gene4':False, 'gene5':False, 'gene12':True, 'gene6':True, 'gene7':True, 'gene8':True, 'gene9':False, 'gene10':False, 'gene11':False}
        gene_sequence = {'gene1':'Sequence1', 'gene2':'Sequence1', 'gene3':'Sequence1', 'gene4':'Sequence1', 'gene5':'Sequence1', 'gene12':'Sequence1', 'gene6':'Sequence2', 'gene7':'Sequence2', 'gene8':'Sequence2', 'gene9':'Sequence2', 'gene10':'Sequence2', 'gene11':'Sequence3'}
        index_file = {'gene1':0, 'gene2':1, 'gene3':2, 'gene4':3, 'gene5':4, 'gene6':6, 'gene7':7,'gene12':5, 'gene8':8, 'gene9':9, 'gene10':10, 'gene11':11}
        
        g = nx.Graph(nx.drawing.nx_pydot.read_dot(input_file))
        G = g.subgraph([gene for gene in g.nodes() if gene in index_file])
        h = nx.Graph(nx.drawing.nx_pydot.read_dot(expected_graph_1))
        self.assertTrue(nx.is_isomorphic(G,h))
        
        closest_genes_dict = {}
        for gene in G.nodes():
            end_list = gg.ends_of_sequence(gene, G, gene_present, gene_sequence)
            for i in end_list:
                closest_genes_dict[i] = gg.closest_gene(i,G, gene_present, gene_sequence)
        self.assertDictEqual(closest_genes_dict, {'gene1':False, 'gene3':False, 'gene7':False, 'gene8':False})

    
    
