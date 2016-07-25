import networkx as nx
from madansi.BlastHit import BlastHit
from madansi.GenePresent import GenePresent
from collections import deque

class WalkGraphs(object):
    
    def __init__(self,graphfile,filteredfile, outputgraphfile):
        self.graphfile = graphfile
        self.filteredfile = filteredfile
        self.outputgraphfile = outputgraphfile
        
    def open_graph_file(self):
        """Open the given graph file for searching"""
        try:
            g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
            return g
        except IOError:
            raise IOError("Error opening this file")
    
    
    def add_node_attribute(self):
        """Adds node attribute to the graph based on whether the gene is given as present in the lookup table as well as the contig that the gene is in"""    
        g = self.open_graph_file()
        gene_dict = GenePresent.construct_dictionary(self)
        for gene in nx.nodes_iter(g):
            g.node[gene]['Contig'] = self.find_sequence(gene)
            if gene_dict[gene]:
                g.node[gene]['present']=True
            else:
                g.node[gene]['present']=False
        return g
    
    def create_subgraph(self):
        """Creates a subgraph with nodes representing the genes that are marked as present"""
        g = self.add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['present']])
        return h
    
    def starting_gene(self):
        """Chooses a starting gene and returns it"""
        h = self.create_subgraph()
        for node in nx.nodes_iter(h):
            start_gene = node
            break
        return start_gene
    
    def find_sequence(self, gene):
        """Given a gene will find out the sequence that it is in"""
        try:
            with open(self.filteredfile,'r') as f:
                for line in f:
                    l = line.rstrip().split('\t')
                    if l[0] == gene:
                        return l[1]
                f.close()
        except IOError:
            raise IOError("Error opening this file")
       
    def construct_sequence_list(self):
        """Constructs a list that will contain all of the names of the sequences"""
        sequence_list =[]
        try:
            with open(self.filteredfile,'r') as f:
                for line in f:
                    l = line.rstrip().split('\t')
                    if l[1] not in sequence_list:
                       sequence_list.append(l[1])
                return sequence_list
                f.close()
                
        except IOError:
            raise IOError("Error opening this file")
        
    def find_ends_of_sequence(self,gene):
        """Given a gene, will find the ends of the sequence containing that gene"""
        g = self.open_graph_file()
        h = self.create_subgraph()
        try:
            end_list = []
            sequence = h.node[gene]['Contig']
            for node in nx.nodes_iter(h):
                if h.node[node]['Contig'] == sequence:
                    neighbor_list=[]
                    for neighbor in h.neighbors(node):
                        if h.node[neighbor]['Contig'] == sequence:
                            neighbor_list.append(neighbor)
                    if len(neighbor_list)==1:
                        end_list.append(node)
            return end_list
        except KeyError:
            raise KeyError('Given gene is not present')
        
    def order_sequences(self,start_gene): #Need to add another method to allow for possible reorientations
        """Constructs a generator to find the closest neighbor from one end""" 
        end_list = self.find_ends_of_sequence(start_gene)
        g = self.add_node_attribute()
   
        neighbors = g.neighbors_iter
        
        if start_gene in end_list:
            gene = start_gene
        else:
            gene = end_list[0]
        
        visited = [gene]
        queue = deque([(gene, neighbors(gene))])
        while queue:
            parent, children = queue[0]
            try:
                child = next(children)
                
                if child not in visited:
                    if g.node[child]['present'] and g.node[child]['Contig'] != g.node[gene]['Contig']:
                        yield child
                        break
                    else:
                        yield parent, child
                        visited.append(child)
                        queue.append((child, neighbors(child)))
            except StopIteration:
                queue.popleft()
    
    def closest_gene(self,start_gene):
        """From the generator defined in order_sequences extracts the closest gene"""
        output_list = self.order_sequences(start_gene)
        x = output_list.__next__()
        while type(x) == tuple:
            x = output_list.__next__()
        return x
        
    def dictionary_pairs_closest_genes(self):
        """Constructs a dictionary to pair all of the closest genes on separate sequences"""
        closest_genes_dict = {}
        
        g = self.add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['present']])
        
        for node in nx.nodes_iter(h):
            end_list = self.find_ends_of_sequence(node)
            closest_genes_dict[end_list[0]]= self.closest_gene(end_list[0])
            closest_genes_dict[end_list[1]]= self.closest_gene(end_list[1])
        
        visited_genes = []
        
        for end_gene in list(closest_genes_dict.keys()):
            if end_gene not in visited_genes:
                sequence_1 = self.find_sequence(end_gene)
                other_end = [i for i in self.find_ends_of_sequence(end_gene) if i!=end_gene][0]
                sequence_2 = self.find_sequence(other_end)
                visited_genes.append(end_gene)
                visited_genes.append(other_end)
                if sequence_1 == sequence_2:
                    if nx.shortest_path_length(g,end_gene,closest_genes_dict[end_gene]) > nx.shortest_path_length(g,other_end, closest_genes_dict[other_end]):
                        closest_genes_dict[end_gene] = None
                    else:
                        closest_genes_dict[other_end] = None

        return closest_genes_dict        
            
    def ordering_sequences(self):
        """Puts the sequences in order and orientates them as necessary"""
        closest_genes_dict = self.dictionary_pairs_closest_genes()
        sequence_list = self.construct_sequence_list()
        
        visited = []
        order_genes_visited = []
        
        for key in list(closest_genes_dict.keys()):
            if closest_genes_dict[key]==None:
                start_end = key
                break
        
        start_sequence = self.find_sequence(start_end)
        visited.append(start_sequence)
        order_genes_visited.append(start_end)
        ends_list = self.find_ends_of_sequence(start_end)
        other_end = [ i for i in ends_list if i!= start_end][0]
        order_genes_visited.append(other_end)
        
        while set(visited) != set(sequence_list):
            if closest_genes_dict[other_end] != None:
                start_end = closest_genes_dict[other_end]
                visited.append(self.find_sequence(start_end))
                order_genes_visited.append(start_end)
                ends_list = self.find_ends_of_sequence(start_end)
                other_end = [i for i in ends_list if i!= start_end][0]
                order_genes_visited.append(other_end)
        
        return order_genes_visited
            
    def create_linear_subgraph(self):
        """Creates a separate linear subgraph that will just consist of the sequences with the shortest path between them"""
        closest_genes_dict = self.dictionary_pairs_closest_genes()
        for key in list(closest_genes_dict.keys()):
            if closest_genes_dict[key]==None:
                start_end = key
                break
                
        g = self.add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['present']])
        
        for key in list(closest_genes_dict.keys()):
            if closest_genes_dict[key]!= None:
                gene_path = nx.shortest_path(g,key, closest_genes_dict[key])
                list_edges = [(gene_path[i], gene_path[i+1]) for i in range(len(gene_path)  - 1)]
                h.add_nodes_from(gene_path)
                h.add_edges_from(list_edges)
        
        nx.Graph(nx.drawing.nx_pydot.write_dot(h,self.outputgraphfile))