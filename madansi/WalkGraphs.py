import networkx as nx
from madansi.BlastHit import BlastHit
from madansi.GenePresent import GenePresent
from collections import deque
from types import GeneratorType

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
        """Adds node attribute to the graph based on whether the gene is given as Present in the lookup table as well as the contig that the gene is in"""    
        g = self.open_graph_file()
        gene_present_dict = GenePresent.construct_dictionary(self)
        for gene in nx.nodes_iter(g):
            g.node[gene]['Sequence'] = self.find_sequence(gene)
            if gene_present_dict[gene]:
                g.node[gene]['Present']=True
            else:
                g.node[gene]['Present']=False
        return g
    
    def create_subgraph(self):
        """Creates a subgraph with nodes rePresenting the genes that are marked as Present"""
        g = self.add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['Present']])
        return h
    
    def starting_gene(self):
        """Chooses a starting gene and returns it"""
        h = self.create_subgraph()
        if len(h.nodes())>0:
            start_gene = h.nodes()[0]
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
        g = self.add_node_attribute()
        
        for node in nx.nodes_iter(g):
            if g.node[node]['Sequence'] not in sequence_list:
                sequence_list.append(g.node[node]['Sequence'])
        
        return sequence_list
        
    def find_ends_of_sequence(self,gene):
        """Given a gene, will find the ends of the sequence containing that gene"""
        g = self.open_graph_file()
        h = self.create_subgraph()
        try:
            end_list = []
            sequence = h.node[gene]['Sequence']
            for node in nx.nodes_iter(h):
                if h.node[node]['Sequence'] == sequence:
                    neighbor_list=[]
                    for neighbor in h.neighbors(node):
                        if h.node[neighbor]['Sequence'] == sequence:
                            neighbor_list.append(neighbor)
                    if len(neighbor_list)==1:
                        end_list.append(node)
            return end_list
        except KeyError:
            raise KeyError('Given gene is not present')
  
        
    def order_sequences(self,start_gene): 
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
                    if g.node[child]['Present'] and g.node[child]['Sequence'] != g.node[gene]['Sequence']:
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
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['Present']])
        
        for node in nx.nodes_iter(h):
            end_list = self.find_ends_of_sequence(node)
            closest_genes_dict[end_list[0]]= self.closest_gene(end_list[0])
            closest_genes_dict[end_list[1]]= self.closest_gene(end_list[1])
        
        visited_genes = []
        list_end_genes = list(closest_genes_dict.keys())
        
        if len(list_end_genes) == 0:
            pass
        elif len(list_end_genes) == 1:
            for gene in list_end_genes:
                closest_genes_dict[gene] = None
        else:
            sequence_list = [self.find_sequence(gene) for gene in ( list_end_genes[0] or list_end_genes[1])]
            if sequence_list[0]==sequence_list[1] and len(list_end_genes) == 2:
                closest_genes_dict[list_end_genes[0]]= None
                closest_genes_dict[list_end_genes[1]]= None
            else:
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
        sequence_list = self.construct_sequence_list()
        unused_sequences=sequence_list
        g = self.add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['Present']])
        
        sequences_genes_present = []
        for node in h.nodes():
            if h.node[node]['Sequence'] not in sequences_genes_present:
                sequences_genes_present.append(h.node[node]['Sequence'])
            if len(sequences_genes_present) > 1:
                break
        if len(sequence_list) == 0: #No sequences present
            unused_sequences = []
            k = nx.Graph()
            nx.drawing.nx_pydot.write_dot(k,self.outputgraphfile)
        elif len(sequence_list) == 1: #Only one sequence present
            nx.drawing.nx_pydot.write_dot(g,self.outputgraphfile)
            unused_sequences = []
        elif len(sequences_genes_present) == 0: #More than one sequence with no genes on them present
            k = nx.Graph()
            nx.drawing.nx_pydot.write_dot(k,self.outputgraphfile)
        elif len(sequences_genes_present) == 1: #Exactly one sequence with at least one gene on
            sequence_present = sequences_genes_present[0]
            unused_sequences.remove(sequence_present)
            k = g.subgraph([gene for gene in g.nodes() if g.node[gene]['Sequence'] == sequence_present])
            nx.drawing.nx_pydot.write_dot(k,self.outputgraphfile)
        else: #At least two sequences with genes on
        
            closest_genes_dict = self.dictionary_pairs_closest_genes()
            
            if len(closest_genes_dict) == 0:
                if len(h.nodes())!=0:
                    example_present_gene = h.nodes()[0]
                    sequence = h.node[example_present_gene]['Sequence']
                    h = g.subgraph([node for node in g.nodes() if g.node[node]['Sequence']==sequence])
                    
                
            else:
                for key in list(closest_genes_dict.keys()):
                    if closest_genes_dict[key]==None:
                        start_end = key
                        break
                   
                unused_sequences=sequence_list
                for node in nx.nodes_iter(h):
                    if h.node[node]['Sequence'] in unused_sequences:
                        unused_sequences.remove(h.node[node]['Sequence'])
                   
                   
                for key in list(closest_genes_dict.keys()):                
                    if closest_genes_dict[key]!= None:
                        gene_path = nx.shortest_path(g,key, closest_genes_dict[key])
                        list_edges = [(gene_path[i], gene_path[i+1]) for i in range(len(gene_path)  - 1)]
                                           
                        for node in gene_path:
                            h.add_node(node)
                            h.node[node]['Present'] = g.node[node]['Present']
                            h.node[node]['Sequence'] = g.node[node]['Sequence']
                            if h.node[node]['Sequence'] in unused_sequences:
                                unused_sequences.remove(h.node[node]['Sequence'])
                                   
                        for edge in list_edges:
                            if 'weight' in g.edge[edge[0]][edge[1]]:
                                w = g.edge[edge[0]][edge[1]]['weight']
                                h.add_edge(edge[0], edge[1],weight=w)
                            else:
                                h.add_edge(edge[0],edge[1])
            nx.drawing.nx_pydot.write_dot(h,self.outputgraphfile)   
        return unused_sequences
        

        
        
        
        