import networkx as nx
from madansi.BlastHit import BlastHit
from madansi.GenePresent import GenePresent
from madansi.DepthFirstSearch import DepthFirstSearch

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
    
    def open_filtered_file(self):
        """Opens the fitlered file for searching"""
        try:
            f = open(self.filteredfile,'r')
            return f
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
        f = self.open_filtered_file()
        for line in f:
            l = line.rstrip().split('\t')
            if l[0] == gene:
                return l[1]
    
    def construct_contig_list(self):
        """Constructs a list that will contain all of the names of the contigs"""
        contig_list =[]
        f = self.open_filtered_file()
        for line in f:
            l = line.rstrip().split('\t')
            if l[1] not in contig_list:
               contig_list.append(l[1])
        return contig_list
        
    def find_neighbors_on_contig(self, gene):
        """Given a gene in the graph checks whether it is marked and if it is returns a list of all its
        neighbors also on the graph"""
        g = self.open_graph_file()
        h = self.create_subgraph()
        if h.node[gene]:
            contig = self.find_contig(gene)
            contig_neighbor_list = []
            for neighbor in g.neighbors(gene):
                if g.node[neighbor]['present'] and self.find_contig(neighbor) == contig:
                    contig_neighbor_list.append(neighbor)
        else:
            raise KeyError('Given gene not marked as present')
        return contig_neighbor_list
    
    def start_and_ends_contigs(self, gene):
        """Finds the genes at either end of the contig, defined from a marked gene on it, where the ends
         are defined to be based upon where the ends of the contig are marked on the graph"""
        contig_neighbor_list = self.find_neighbors_on_contig()
     
    def ordering_contigs(self):
        """Puts the contigs in order and orientates them as necessary"""
        visited = []
        contig_list = self.construct_contig_list()
        start_gene = self.starting_gene()
        start_contig = self.find_contig(start_gene)
        visited.append(start_contig)
        while set(visited) != set(contig_list):
            
        
        
#First find a contig that is marked as present in the graph file.
#Then look at all nodes that are neighbors of this node and see if there is a contig in the same or a different sequence (should be marked as being present in the lookup table) here. If it is in the same sequence - which we will know from looking at the filtered fasta file and the contig name- so should be the second column given that we have swapped the columns initially.
#If a second contig has not been found then we would need to extend our search to two neighbours away and continue in that way. However, if a contig has been found then we need to go to the last present contig in that sequence and continue the search from there.      
   
#Add node attributes to make sure that the genes are actually present.        

#Create subgraph with all the genes present in the contig

#Work with this graph to find an initial gene to start with

#Make a note of the contig that it is from- need to relate this back to the filtered file and should take the second column of the row where this gene is first. 
#Open the filtered file for reading:

#Iterate through the filtered file to find the record that includes the relevant gene and return the second entry of this row

#Starting with the chosen gene will want to look through the neighbours to see if we can find a gene that is marked and belongs to a different contig.

    