import networkx as nx
from madansi.BlastHit import BlastHit
from madansi.GenePresent import GenePresent
from madansi.DepthFirstSearch import DepthFirstSearch

#First find a contig that is marked as present in the graph file.
#Then look at all nodes that are neighbors of this node and see if there is a contig in the same or a different sequence (should be marked as being present in the lookup table) here. If it is in the same sequence - which we will know from looking at the filtered fasta file and the sample name- so should be the second column given that we have swapped the columns initially.
#If a second sample has not been found then we would need to extend our search to two neighbours away and continue in that way. However, if a sample has been found then we need to go to the last present contig in that sequence and continue the search from there.

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
            raise Error("Error opening this file")

    
# Choose starting sequence and will need to iterate through all of the sequences until finished.
    sequence_list=[]
    visited=[]       
         
        
    def create_subgraph(self):
        """Creates a subgraph with nodes representing the genes that are marked as present"""
        g = DepthFirstSearch(self.graphfile,self.filteredfile).add_node_attribute()
        h = g.subgraph([gene for gene in g.nodes() if g.node[gene]['present']])
        return h
    
    def starting_gene(self):
        """Chooses a starting gene and returns it"""
        h = self.create_subgraph()
        start_gene = h.node[0]
        return start_gene
    
    def find_sequence(self, gene):
        """Given a gene will find out the sequence that it is in"""
        f = open(self.filteredfile,'r')
        for line in f:
            l = line.rstrip().split('\t')
            if l[0] == gene:
                return l[1]
         
         
#Add node attributes to make sure that the genes are actually present.        

#Create subgraph with all the genes present in the sample

#Work with this graph to find an initial gene to start with

#Make a note of the sample that it is from- need to relate this back to the filtered file and should take the second column of the row where this gene is first. 
#Open the filtered file for reading:

#Iterate through the filtered file to find the record that includes the relevant gene and return the second entry of this row

    