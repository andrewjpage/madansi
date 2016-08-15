from madansi.ContigSearching                    import ContigSearching
from madansi.RefineContigNeighbours             import RefineContigNeighbours
from madansi.ContigGraph                        import ContigGraph
from madansi.IterateJoiningContigComponents     import IterateJoiningContigComponents

class ProduceOrderedContigGraph(object):
    
    def __init__(self, gene_detector, filtered_graph, filtered_blast_hits_file, output_refined_contig_graph):
        self.gene_detector                  = gene_detector
        self.filtered_graph                 = filtered_graph
        self.filtered_blast_hits_file       = filtered_blast_hits_file
        self.output_refined_contig_graph    = output_refined_contig_graph
    
    def produce_ordered_contig_graph(self):
        contig_searching = ContigSearching(self.gene_detector, self.filtered_graph)
        contig_searching.expand_all_contigs()
        
        refine_neighbouring_contigs = RefineContigNeighbours(contig_searching.neighbouring_contigs, self.filtered_graph, self.filtered_blast_hits_file, self.gene_detector)
        refine_neighbouring_contigs.refine_contig_neighbours()
        
        contig_graph_refined    = ContigGraph(refine_neighbouring_contigs.refined_neighbouring_contigs)
        contig_graph_unrefined  = ContigGraph(contig_searching.neighbouring_contigs)
        graph_refined           = contig_graph_refined.create_contig_subgraph()
        graph_unrefined         = contig_graph_unrefined.create_contig_subgraph()
        
        iterate_joining_contig_components   = IterateJoiningContigComponents(graph_unrefined, self.output_refined_contig_graph)
        ordered_contig_graph                = iterate_joining_contig_components.iterate_joining(graph_refined)
        iterate_joining_contig_components.output_graph(ordered_contig_graph)
        return ordered_contig_graph
