#!/usr/bin/env python3

import argparse
import pprint
import os
from madansi.ContigSearching import ContigSearching
from madansi.GeneDetector import GeneDetector
from madansi.GraphParser import GraphParser
from madansi.ContigGraph import ContigGraph
from madansi.FilterBlastComparison import FilterBlastComparison
from madansi.ContigSpanningTree import ContigSpanningTree
from madansi.UnusedContigs import UnusedContigs
from madansi.RefineContigNeighbours import RefineContigNeighbours
from madansi.RemoveHeaviestEdges import RemoveHeaviestEdges

parser = argparse.ArgumentParser(description = 'Script to take in the input assembly file, blast hits file and graph file and output a list of how each of the \
contigs are connected to each other and the number of iterations it takes to get between the two.')
parser.add_argument('input_assembly_file', help = 'Path to the input assembly file.', type = str)
parser.add_argument('blast_hits_file', help = 'Path to the blast hits file.', type = str)
parser.add_argument('filtered_blast_hits_file', help = 'Path to the output filtered blast hits file.', type = str)
parser.add_argument('graph_file', help = 'Path to the inout graph file.', type = str)
parser.add_argument('output_contig_graph_file', help = 'Path to the output contig graph file. This should be a dot file.', type = str)
parser.add_argument('output_spanning_tree', help = 'Path to the output spanning tree of contigs', type = str)
parser.add_argument('output_contig_file', help = 'Path to the file containing the names of the unused contigs', type = str)
args = parser.parse_args()



filtered_blast_hits_file = FilterBlastComparison(args.blast_hits_file, args.filtered_blast_hits_file,  bit_score=200)
filtered_blast_hits_file.filter()

graph_parser = GraphParser(args.graph_file)
filtered_graph = graph_parser.graph

gene_detector = GeneDetector(args.input_assembly_file, args.filtered_blast_hits_file)
gene_detector.contigs_to_genes()

unused_contigs = UnusedContigs(gene_detector, gene_detector.assembly.sequence_names(), args.output_contig_file )
unused_contigs.contigs_not_in_filtered_file()
unused_contigs.output_unused_contigs()

contig_searching = ContigSearching(gene_detector, filtered_graph)
contig_searching.expand_all_contigs()

refine_neighbouring_contigs = RefineContigNeighbours(contig_searching.neighbouring_contigs, filtered_graph, args.filtered_blast_hits_file)
refine_neighbouring_contigs.refine_contig_neighbours()

contig_graph = ContigGraph(refine_neighbouring_contigs.refined_neighbouring_contigs, args.output_contig_graph_file)
contig_graph.create_contig_subgraph()
contig_graph.output_contig_graph()

#contig_graph = ContigGraph(refine_neighbouring_contigs.refine_contig_neighbours(), 'temp_contig_graph')


#rhe = RemoveHeaviestEdges(contig_graph.create_contig_subgraph(), args.output_contig_graph_file)
#rhe.output_simplified_graph()


#contig_spanning_tree = ContigSpanningTree(contig_graph.contig_graph, args.output_spanning_tree)
#contig_spanning_tree.construct_spanning_tree()
 