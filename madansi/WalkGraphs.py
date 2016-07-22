import networkx as nx
from madansi.BlastHit import BlastHit
from madansi.GenePresent import GenePresent
from madansi.DepthFirstSearch import DepthFirstSearch
from collections import deque

class WalkGraphs(object):
    
    def __init__(self,graphfile,filteredfile):
        self.graphfile = graphfile
        self.filteredfile = filteredfile
        
    def open_graph_file(self):
        """Open the given graph file for searching"""
        try:
            g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
            return g
        except IOError:
            raise IOError("Error opening this file")
    
    def create_subgraph(self):
        """Creates a subgraph with nodes representing the genes that are marked as present"""
        g = DepthFirstSearch(self.graphfile,self.filteredfile).add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['present']])
        return h
    
    def starting_gene(self):
        """Chooses a starting gene and returns it"""
        h = self.create_subgraph()
        for node in nx.nodes_iter(h):
            start_gene = node
            break
        return start_gene
    
    def find_contig(self, gene):
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
       
    def construct_contig_list(self):
        """Constructs a list that will contain all of the names of the contigs"""
        contig_list =[]
        try:
            with open(self.filteredfile,'r') as f:
                for line in f:
                    l = line.rstrip().split('\t')
                    if l[1] not in contig_list:
                       contig_list.append(l[1])
                return contig_list
                f.close()
                
        except IOError:
            raise IOError("Error opening this file")
        
    def find_ends_of_contig(self,gene):
        """Given a gene, will find the ends of the contig containing that gene"""
        g = self.open_graph_file()
        h = self.create_subgraph()
        try:
            end_list = []
            contig = h.node[gene]['Contig']
            for node in nx.nodes_iter(h):
                if h.node[node]['Contig'] == contig:
                    neighbor_list=[]
                    for neighbor in h.neighbors(node):
                        if h.node[neighbor]['Contig'] == contig:
                            neighbor_list.append(neighbor)
                    if len(neighbor_list)==1:
                        end_list.append(node)
            return end_list
        except KeyError:
            raise KeyError('Given gene is not present')
        
    def order_contigs(self,start_gene): #Need to add another method to allow for possible reorientations
        """Constructs a generator to find the closest neighbor from one end""" 
        end_list = self.find_ends_of_contig(start_gene)
        g = DepthFirstSearch(self.graphfile,self.filteredfile).add_node_attribute()
   
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
        """From the generator defined in order_contigs extracts the closest gene"""
        output_list = self.order_contigs(start_gene)
        x = output_list.__next__()
        while type(x) == tuple:
            x = output_list.__next__()
        return x
        
        
        
    def dictionary_pairs_closest_genes(self):
        """Constructs a dictionary to pair all of the closest genes on separate contigs"""
        end_genes_dict = {}
        
        g = DepthFirstSearch(self.graphfile,self.filteredfile).add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['present']])
        
        for node in nx.nodes_iter(h):
            end_list = self.find_ends_of_contig(node)
            #Maybe want to change this so that it doesn't define stuff twice 
            end_genes_dict[end_list[0]]= self.closest_gene(end_list[0])
            end_genes_dict[end_list[1]]= self.closest_gene(end_list[1])
        
        return end_genes_dict
        
            
  #  def ordering_contigs(self):
  #      """Puts the contigs in order and orientates them as necessary"""
  #      visited = []
  #      contig_list = self.construct_contig_list()
  #      start_gene = self.starting_gene()
  #      start_contig = self.find_contig(start_gene)
  #      visited.append(start_contig)
  #      while set(visited) != set(contig_list):
            
        
    