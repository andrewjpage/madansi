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
        
    def order_contigs(self):
        start_gene = self.starting_gene()
        end_list = self.find_ends_of_contig(start_gene)
        g = DepthFirstSearch(self.graphfile,self.filteredfile).add_node_attribute()
                
        neighbors = g.neighbors_iter
        
        visited = set(end_list[0])
        queue = deque([(end_list[0], neighbors(end_list[0]))])
        while queue:
            parent, children = queue[0]
            try:
                child = next(children)
                if child not in visited:
                    if not (g.node[child]['present'] and g.node[child]['Contig'] != contig):
                        yield parent, child
                        visited.add(child)
                        queue.append((child, neighbors(child)))
                    else:
                        yield child
            except StopIteration:
                queue.popleft()
        
     
  #  def ordering_contigs(self):
  #      """Puts the contigs in order and orientates them as necessary"""
  #      visited = []
  #      contig_list = self.construct_contig_list()
  #      start_gene = self.starting_gene()
  #      start_contig = self.find_contig(start_gene)
  #      visited.append(start_contig)
  #      while set(visited) != set(contig_list):
            
        
        
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

    