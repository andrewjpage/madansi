import networkx as nx

class ContigSubgraph(object):
    def __init__(self, graph, output_graph_file):
        self.graph = graph
        self.output_graph_file = output_graph_file
    
    def create_subgraph(self):
        nodes_degree_at_most_two = []
        for node in self.graph.nodes():
            if self.graph.degree(node)<=2:
                nodes_degree_at_most_two.append(node)

        edge_list = self.graph.edges(nodes_degree_at_most_two, data=True)
        contig_subgraph = nx.Graph()
        contig_subgraph.add_edges_from(edge_list)
        return contig_subgraph
    
    def output_subgraph(self):
        contig_subgraph = self.create_subgraph()
        nx.drawing.nx_pydot.write_dot(contig_subgraph, self.output_graph_file)