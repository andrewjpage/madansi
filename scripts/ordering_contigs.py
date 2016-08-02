#!/usr/bin/env python3

import argparse
import pprint
from madansi.ContigSearching import ContigSearching
from madansi.GeneDetector import GeneDetector
from madansi.GraphParser import GraphParser
from madansi.ContigGraph import ContigGraph
from madansi.FilterBlastComparison import FilterBlastComparison

parser = argparse.ArgumentParser(description = 'Script to take in the input assembly file, blast hits file and graph file and output a list of how each of the \
contigs are connected to each other and the number of iterations it takes to get between the two.')
parser.add_argument('input_assembly_file', help = 'Path to the input assembly file.', type = str)
parser.add_argument('blast_hits_file', help = 'Path to the blast hits file.', type = str)
parser.add_argument('filtered_blast_hits_file', help = 'Path to the output filtered blast hits file.', type = str)
parser.add_argument('graph_file', help = 'Path to the inout graph file.', type = str)
parser.add_argument('output_contig_graph_file', help = 'Path to the output contig graph file. This should be a dot file.', type = str)
args = parser.parse_args()


graph_parser = GraphParser(args.graph_file)
filtered_graph = graph_parser.graph

filtered_blast_hits_file = FilterBlastComparison(args.blast_hits_file, args.filtered_blast_hits_file,  bitscore=200)
filtered_blast_hits_file.filter()

gene_detector = GeneDetector(args.input_assembly_file, args.filtered_blast_hits_file)
gene_detector.contigs_to_genes()

contig_searching = ContigSearching(gene_detector, filtered_graph)
contig_searching.expand_all_contigs()

contig_graph = ContigGraph(contig_searching.neighbouring_contigs, args.output_contig_graph_file)
contig_graph.create_contig_subgraph()

