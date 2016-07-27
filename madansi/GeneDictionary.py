import networkx as nx

class GeneDictionary(object):
    
    def __init__(self,graphfile, filteredfile):
        self.graphfile = graphfile
        self.filteredfile = filteredfile
    
    def open_graph_file(self):
        """Open the given graph file for searching"""
        try:
            g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
            return g
        except IOError:
            raise IOError("Error opening this file")
    

    def index_filtered_file(self):
        """Constructs a dictionary linking the name of the gene with the line of the filtered file that it is on"""
        index_lines = {}
        count = 0
    
        try:
            f = open(self.filteredfile)
        except IOError:
            raise Error("Error opening this file")
        
        for line in f:
            l = line.split("\t")
            index_lines[l[0]] = count
            count += 1
    
        f.close()
        return index_lines
        
    def construct_dictionary(self):
        gene_dict = {}
        g = self.open_graph_file()
        
        
            
        for node in g.nodes():
            try:
                line_blast_file = self.index_filtered_file()[node]
            except KeyError:
                raise KeyError("Gene not present in this file")
            
            try:
                f = open(self.filteredfile)
            except IOError:
                raise Error("Error opening this file")
                
            for i,line in enumerate(f):
                if i == line_blast_file:
                    
                    l = line.split("\t")
                    if l[8] <= l[9]:
                        orientation = True
                    else:
                        orientation = False
                        
                    sequence = l[1]
                    
                    if float(l[11])>= 150:
                        present = True
                    else:
                        present = False
            f.close()
            gene_dict[node]= (present,sequence,orientation)

        return gene_dict
            