import networkx as nx

class ContigOrientation(object):
    def __init__(self,contig_graph, gene_detector):
        self.contig_graph   = contig_graph
        self.gene_detector  = gene_detector
        self.contigs = self.gene_detector.contigs
    
    def find_contig_orientation(self, contig):
        forward = 0
        reverse = 0
        for gene in self.contigs[contig].gene_objects.keys():
            if self.contigs[contig].gene_objects[gene].orientation > 0:
                forward += 1
            else:
                reverse += 1
        if forward >= reverse:
            return 1
        else:
            return -1
        
    def orientate_edges(self, contig, contig_list):
        visited_contigs = [contig]
        queue = [contig]
        while sorted(visited_contigs) != sorted(contig_list):
            print(visited_contigs)
            print(queue)
            contig = queue[0]
            for neighbour in self.contig_graph.neighbors(contig):
                if neighbour not in visited_contigs and neighbour not in queue:
                    queue.append(neighbour)
            for edge in self.contig_graph.edges(contig):
                contig_1_orientation = self.find_contig_orientation(edge[0])
                contig_2_orientation = self.find_contig_orientation(edge[1])
                if contig_1_orientation == 1 and contig_2_orientation == 1: #forward - forward
                    self.contig_graph.add_edge(edge[0], edge[1],weight= 1)
                elif contig_1_orientation == 1 and contig_2_orientation == -1: #forward - reverse
                    self.contig_graph.add_edge(edge[0], edge[1],weight=2)
                elif contig_1_orientation == -1 and contig_2_orientation == 1: #reverse - forward
                    self.contig_graph.add_edge(edge[0], edge[1],weight = 3)
                else:
                    self.contig_graph.add_edge(edge[0], edge[1],weight = 4) #reverse - reverse
                if edge[1] not in visited_contigs:
                    visited_contigs.append(edge[1])
            queue.remove(contig)

    
    def repeat_all_connected_components(self):
        for connected_component in sorted(nx.connected_components(self.contig_graph)):
            start_contig = sorted(connected_component)[0]
            self.orientate_edges(start_contig, sorted(connected_component))
        return self.contig_graph
    
    def output_graph(self, output_graph_file):
        nx.drawing.nx_pydot.write_dot(self.contig_graph, output_graph_file)