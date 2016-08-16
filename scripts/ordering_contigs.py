#!/usr/bin/env python3

import argparse
import os
from madansi.ContigSearching                import ContigSearching
from madansi.GeneDetector                   import GeneDetector
from madansi.GraphParser                    import GraphParser
#from madansi.ContigGraph                    import ContigGraph
from madansi.FilterBlastComparison          import FilterBlastComparison
from madansi.UnusedContigs                  import UnusedContigs
#from madansi.RefineContigNeighbours         import RefineContigNeighbours
#from madansi.JoiningContigComponents        import JoiningContigComponents
#from madansi.IterateJoiningContigComponents import IterateJoiningContigComponents 
#from madansi.ContigOrientation              import ContigOrientation
from madansi.ProduceOrderedContigGraph      import ProduceOrderedContigGraph
from madansi.GraphToFasta                   import GraphToFasta

parser = argparse.ArgumentParser(description =  'Script to take in the input assembly file, blast hits file and graph file and output a list of how each of the \
                                                contigs are connected to each other and the number of iterations it takes to get between the two.')
parser.add_argument('input_assembly_file',          help = 'Path to the input assembly file.',                                          type = str)
parser.add_argument('blast_hits_file',              help = 'Path to the blast hits file.',                                              type = str)
parser.add_argument('filtered_blast_hits_file',     help = 'Path to the output filtered blast hits file.',                              type = str)
parser.add_argument('input_graph_file',             help = 'Path to the input graph file.',                                             type = str)
parser.add_argument('output_refined_contig_graph',  help = 'Path to the output refined contig graph file. This should be a dot file.',  type = str)
parser.add_argument('output_contig_file',           help = 'Path to the file containing the names of the unused contigs.',              type = str)
parser.add_argument('output_fasta_file',                                                                                                type = str)
args = parser.parse_args()

filtered_blast_hits_file = FilterBlastComparison(args.blast_hits_file, args.filtered_blast_hits_file,  bit_score=200)
filtered_blast_hits_file.filter()

graph_parser    = GraphParser(args.input_graph_file)
filtered_graph  = graph_parser.graph

gene_detector   = GeneDetector(args.input_assembly_file, args.filtered_blast_hits_file)
gene_detector.contigs_to_genes()

unused_contigs  = UnusedContigs(gene_detector, gene_detector.assembly.sequence_names(), args.output_contig_file )
unused_contigs.contigs_not_in_filtered_file()

produced_ordered_contig_graph   = ProduceOrderedContigGraph(gene_detector, filtered_graph, args.filtered_blast_hits_file, args.output_refined_contig_graph)
ordered_contig_graph            = produced_ordered_contig_graph.produce_ordered_contig_graph()
contig_ends                     = produced_ordered_contig_graph.contig_ends

graph_to_fasta = GraphToFasta(args.input_assembly_file, ordered_contig_graph, args.output_fasta_file, contig_ends)
graph_to_fasta.create_fasta_file()

#contig_searching = ContigSearching(gene_detector, filtered_graph)
#contig_searching.expand_all_contigs()
#
#refine_neighbouring_contigs = RefineContigNeighbours(contig_searching.neighbouring_contigs, filtered_graph, args.filtered_blast_hits_file, gene_detector)
#refine_neighbouring_contigs.refine_contig_neighbours()
#
#contig_graph_refined    = ContigGraph(refine_neighbouring_contigs.refined_neighbouring_contigs)
#contig_graph_unrefined  = ContigGraph(contig_searching.neighbouring_contigs)
#graph_refined           = contig_graph_refined.create_contig_subgraph()
#graph_unrefined         = contig_graph_unrefined.create_contig_subgraph()
#
#iterate_joining_contig_components   = IterateJoiningContigComponents(graph_unrefined, args.output_refined_contig_graph)
#ordered_contig_graph                = iterate_joining_contig_components.iterate_joining(graph_refined)
#iterate_joining_contig_components.output_graph(ordered_contig_graph)

#unused_contigs.contigs_not_in_refined_graph(ordered_contig_graph)
#unused_contigs.output_unused_contigs()
#
#contig_orientation = ContigOrientation(graph_unrefined, gene_detector)
#contig_orientation.repeat_all_connected_components()
#contig_orientation.output_graph(args.output_refined_contig_graph)
#
