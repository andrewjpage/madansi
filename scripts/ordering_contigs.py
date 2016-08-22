#!/usr/bin/env python3

import argparse
import os
from madansi.RunBLAST                       import RunBLAST
from madansi.ContigSearching                import ContigSearching
from madansi.GeneDetector                   import GeneDetector
from madansi.GraphParser                    import GraphParser
from madansi.FilterBlastComparison          import FilterBlastComparison
from madansi.UnusedContigs                  import UnusedContigs
from madansi.ProduceOrderedContigGraph      import ProduceOrderedContigGraph
from madansi.IterateJoiningContigComponents import IterateJoiningContigComponents
from madansi.GraphToFasta                   import GraphToFasta
from madansi.RemoveContigs                  import RemoveContigs   

parser = argparse.ArgumentParser(description =  'Script to take in the input assembly file, blast hits file and graph file and output a list of how each of the \
                                                contigs are connected to each other and the number of iterations it takes to get between the two.')
parser.add_argument('input_reference',              help ='path to the input reference fasta file',                                     type = str)
parser.add_argument('output_reference',             help = 'path to the output reference fasta file with switched columns',             type = str)
parser.add_argument('output_database',              help = 'path to the output database',                                               type = str)
parser.add_argument('input_assembly_file',          help = 'Path to the input assembly file.',                                          type = str)
parser.add_argument('blast_hits_file',              help = 'Path to the blast hits file.',                                              type = str)
parser.add_argument('filtered_blast_hits_file',     help = 'Path to the output filtered blast hits file.',                              type = str)
parser.add_argument('input_graph_file',             help = 'Path to the input graph file.',                                             type = str)
parser.add_argument('output_refined_contig_graph',  help = 'Path to the output refined contig graph file. This should be a dot file.',  type = str)
parser.add_argument('output_fasta_file',            help = 'Path to the output fasta file with contigs grouped and orientated',         type = str)
args = parser.parse_args()

rb = RunBLAST(args.input_assembly_file, args.input_reference, args.output_reference, args.output_database, args.blast_hits_file)
rb.run_switch_columns_database()
rb.make_reference_database()
rb.run_BLAST()

filtered_blast_hits_file = FilterBlastComparison(args.blast_hits_file, args.filtered_blast_hits_file,  bit_score=200)
filtered_blast_hits_file.filter()

graph_parser    = GraphParser(args.input_graph_file)
filtered_graph  = graph_parser.graph

gene_detector   = GeneDetector(args.input_assembly_file, args.filtered_blast_hits_file)
gene_detector.contigs_to_genes()
gene_detector.assembly.sequence_names()
sequences = gene_detector.assembly.sequences

unused_contigs  = UnusedContigs(gene_detector, args.output_fasta_file, args.input_assembly_file )
unused_contigs.contigs_not_in_filtered_file()

produced_ordered_contig_graph   = ProduceOrderedContigGraph(gene_detector, filtered_graph, args.filtered_blast_hits_file, args.output_refined_contig_graph, sequences)
ordered_contig_graph            = produced_ordered_contig_graph.produce_ordered_contig_graph()
contig_ends                     = produced_ordered_contig_graph.contig_ends

remove_contigs = RemoveContigs(ordered_contig_graph)
ordered_contig_graph_filtered = remove_contigs.remove_extra_contigs()
output_filtered_graph = IterateJoiningContigComponents(ordered_contig_graph_filtered, args.output_refined_contig_graph)
output_filtered_graph.output_graph(ordered_contig_graph_filtered)

graph_to_fasta = GraphToFasta(sequences, ordered_contig_graph_filtered, args.output_fasta_file, contig_ends)
graph_to_fasta.create_fasta_file_combined_contigs()

unused_contigs.contigs_not_in_filtered_graph(ordered_contig_graph_filtered)
unused_contigs.add_unused_contigs_to_end()
