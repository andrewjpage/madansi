import networkx as nx
from madansi.GenePresent import GenePresent
from collections import deque
import sys

class GenerateGraph(object):
	
    def __init__(self,graphfile, filteredfile, outputgraphfile, unused_sequence_file):
        self.graphfile = graphfile
        self.filteredfile = filteredfile
        self.outputgraphfile = outputgraphfile
        self.unused_sequence_file = unused_sequence_file
	
    def open_graph_file(self):
        try:
            g=nx.Graph(nx.drawing.nx_pydot.read_dot(self.graphfile))
            print('Opened graph file')
            return g
        except IOError:
            raise IOError("Error opening this file")
    
    def create_index_file(self):
        index_file = GenePresent.index_filtered_file(self)
        print('Created index of files')
        return index_file
        
    def create_gene_present_dict(self):
        gene_present = GenePresent.construct_dictionary_present(self)
        print('Created dictionary of genes present')
        return gene_present
    
    def create_gene_sequence_dict(self):
        gene_sequence = GenePresent.construct_dictionary_sequence(self)
        print('Created dictionary of genes sequences')
        return gene_sequence
    
    def create_gene_orientation_dict(self):
        gene_orientation = GenePresent.construct_dictionary_orientation(self)
        print('Created dictionary of genes orientation')
        return gene_orientation        
    
    def ends_of_sequence(self,start_gene, graph, gene_present, gene_sequence): #NB potential problem if the sample with genes on is not linear - i.e if it 
        try:
            end_list = []
            sequence = gene_sequence[start_gene]
            visited = []
            original_start_gene = start_gene
            
            if all(gene_sequence[neighbor] != sequence for neighbor in graph.neighbors(start_gene) if gene_present[neighbor]): 
                end_list.append(start_gene)
            
            else: 
                
                while len(end_list)<2:
                        if len([gene for gene in graph.neighbors(start_gene) if gene_sequence[gene] == sequence and gene_present[gene]]) == 1:
                            
                            if start_gene not in end_list:
                                end_list.append(start_gene)
                                visited.append(start_gene)
                                if start_gene != original_start_gene:
                                    start_gene = original_start_gene
                                else:
                                    start_gene = [gene for gene in graph.neighbors(start_gene) if gene_sequence[gene] == sequence and gene_present[gene]][0]
                        else:
                            for i in [gene for gene in graph.neighbors(start_gene) if gene_sequence[gene] == sequence and gene_present[gene]]:

                                if i not in visited:
                                    visited.append(start_gene)
                                    start_gene = i
                                           
                                    break

            return end_list
        except KeyError:
            raise KeyError('Given gene is not present')
            
    def order_sequences(self,start_gene,graph, gene_present, gene_sequence):
        end_list = self.ends_of_sequence(start_gene,graph, gene_present, gene_sequence)

        neighbors = graph.neighbors_iter
        
        if start_gene in end_list:
            gene = start_gene
        else:
            gene = end_list[0]
        
        visited = [gene]
        queue = deque([(gene,neighbors(gene))])
        while queue:
            parent,children = queue[0]
            try:
                child = next(children)
                if child not in visited:
                    if child in gene_present:
                        if gene_present[child] and gene_sequence[child]!=gene_sequence[gene]:
                            yield child
                            break
                        else:
                            yield parent, child
                            visited.append(child)
                            queue.append((child, neighbors(child)))
            except StopIteration:
                queue.popleft()
                

    def closest_gene(self,start_gene, graph, gene_present, gene_sequence):
        output_list = self.order_sequences(start_gene, graph, gene_present, gene_sequence)
        x = output_list.__next__()
        try:
            while type(x) == tuple:
                x = output_list.__next__()
            return x
        except StopIteration:
            return False
        print('Found closest gene')

    
    def generate_graph(self):
        
    #Open the graph file
        g = self.open_graph_file()
            
    #Construct the relevant dictionaries
        index_file = self.create_index_file()
        gene_present = self.create_gene_present_dict()
        gene_sequence = self.create_gene_sequence_dict()
        gene_orientation = self.create_gene_orientation_dict()
                
    #Construct the subgraph of all genes that are contained in the file
        G = g.subgraph([gene for gene in g.nodes() if gene in index_file])
        print('Created subgraph of all genes in the graph that correspond to those in the file')
    
    #Construct the list of sequences present in the file
        sequence_list = []
        for gene in nx.nodes_iter(G):
            if gene in gene_sequence:
                sequence = gene_sequence[gene]
                if sequence not in sequence_list:
                    sequence_list.append(sequence)
        print('Created list of all sequences present')
    
    #Construct list of unused sequences
        unused_sequences = sequence_list
        print('Constructed list of unused sequences')
    
    #Construct subgraph of all the genes that are present
        h = G.subgraph([gene for gene in gene_present if gene_present[gene]])
        print('Constructed subgraph')        
    
    #Construct list of genes present in the file- only important if there are zero or one entries
        sequences_genes_present = []
        for gene in nx.nodes_iter(G):
            if gene in gene_present:
                if gene_present[gene]:
                    if gene_sequence[gene] not in sequences_genes_present:
                        sequences_genes_present.append(gene_sequence[gene])
                    if len(sequences_genes_present) > 1:
                        break
        print('Checked if there were less than 2 sequences present')
    
    #Tests
        
        if len(sequence_list) == 0:
            unused_sequences = []
            k = nx.Graph()
            nx.drawing.nx_pydot.write_dot(k,self.outputgraphfile)
            print('No sequences present')
        elif len(sequence_list) == 1:
            unused_sequences = []
            nx.drawing.nx_pydot.write_dot(G,self.outputgraphfile)
            print('One sequence present')
        elif len(sequences_genes_present) == 0:
            k = nx.Graph()
            nx.drawing.nx_pydot.write_dot(k, self.outputgraphfile)
            print('No genes present')
        elif len(sequences_genes_present) == 1:
            sequence_present = sequences_genes_present[0]
            unused_sequences.remove(sequence_present)
            k = G.subgraph([gene for gene in gene_present if gene_sequence[gene]==sequence_present])
            nx.drawing.nx_pydot.write_dot(k,self.outputgraphfile)
            print('One sequence with genes present')
        else:
            
            closest_genes_dict = self.construct_closest_pairs_dictionary(G, gene_present, gene_sequence)
            
            if len(closest_genes_dict) == 0: #Understand why it is possible to have no elements in closest_genes_dict
                if len(G.nodes()) != 0:
                    example_present_gene = G.nodes()[0]
                    sequence = gene_sequence[example_present_gene]
                    k = G.subgraph([gene for gene in index_file if gene_sequence[gene]==sequence])
                    nx.drawing.nx_pydot.write_dot(k,self.outputgraphfile)
            else:
                for key in list(closest_genes_dict.keys()):
                    if closest_genes_dict[key] == None:
                        start_end = key
                        break
                print('Found starting point')
                
                unused_sequences = sequence_list
                for gene in gene_present:
                    if gene_present[gene]:
                        if gene_sequence[gene] in unused_sequences:
                            unused_sequences.remove(gene_sequence[gene])
                print('Found unused sequences')
               
                for key in list(closest_genes_dict.keys()):
                    if closest_genes_dict[key]!= None:
                                                    
                        if closest_genes_dict[key]==False:
                            end_list = self.ends_of_sequence(key,G, gene_present, gene_sequence)
                            gene_path = nx.shortest_path(G, end_list[0], end_list[1] )
                            for node in gene_path:
                                if node in h.nodes():
                                    h.remove_node(node)
                            if gene_sequence[key] not in unused_sequences:
                                unused_sequences.append(gene_sequence[key])
                            
                        else:
                            gene_path = nx.shortest_path(G, key, closest_genes_dict[key])
                            list_edges = [(gene_path[i], gene_path[i+1]) for i in range(len(gene_path) - 1)]
                            
                            for node in gene_path:
                                h.add_node(node)
                                if gene_sequence[node] in unused_sequences:
                                    unused_sequences.remove(gene_sequence[node])
                            
                            for edge in list_edges:
                                if 'weight' in G.edge[edge[0]][edge[1]]:
                                    w = G.edge[edge[0]][edge[1]]['weight']
                                    h.add_edge(edge[0], edge[1], weight=w)
                                else:
                                    h.add_edge(edge[0],edge[1])
                            
            nx.drawing.nx_pydot.write_dot(h, self.outputgraphfile)
            print('Graph outputted')
        
        with open(self.unused_sequence_file, 'w') as handle:
            for sequence in unused_sequences:
                handle.write(sequence + '\n') 
        handle.close()
        
        return unused_sequences
        
    def construct_closest_pairs_dictionary(self,G, gene_present, gene_sequence):
        closest_genes_dict = {}
        for gene in nx.nodes_iter(G):
            if gene in gene_present:
                if gene_present[gene]:
                    end_list = self.ends_of_sequence(gene, G,gene_present, gene_sequence)
                    for i in end_list:
                        closest_genes_dict[i] = self.closest_gene(i,G, gene_present, gene_sequence)
        
        visited_genes = []
        list_end_genes = list(closest_genes_dict.keys())
        
        if len(list_end_genes) == 0:
            pass
        elif len(list_end_genes) == 1:
            for gene in list_end_genes:
                closest_genes_dict[gene] = None
        else:
            if len(list_end_genes) == 2 and gene_sequence[list_end_genes[0]]==gene_sequence[list_end_genes[1]]:
                closest_genes_dict[list_end_genes[0]] = None
                closest_genes_dict[list_end_genes[1]] = None
            else:
                for end_gene in list(closest_genes_dict.keys()):
                    if end_gene not in visited_genes:
                        sequence_1 = gene_sequence[end_gene]
                        other_end = [i for i in self.ends_of_sequence(end_gene, G, gene_present, gene_sequence) if i!= end_gene][0]
                        sequence_2 = gene_sequence[other_end]
                        visited_genes.append(end_gene)
                        visited_genes.append(other_end)
                        if sequence_1 == sequence_2: #Something to do with a cycle that just considers two sequences that form a cycle.
                            if closest_genes_dict[end_gene] != False and closest_genes_dict[other_end]!= False:
                                if nx.shortest_path_length(G, end_gene, closest_genes_dict[end_gene]) > nx.shortest_path_length(G, other_end, closest_genes_dict[other_end]):
                                    closest_genes_dict[end_gene] = None
                                else:
                                    closest_genes_dict[other_end] = None
        print('Constructed closest_genes_dict')
        
        return closest_genes_dict
                        