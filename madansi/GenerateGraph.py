import networkx as nx
from madansi.GenePresent import GenePresent

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
    
    def ends_of_sequence(self,start_gene, graph, gene_present_dict, gene_sequence):
        try:
            end_list = []
            sequence = gene_sequence[start_gene]
            for node in gene_present_dict:
                if gene_present_dict[node]:
                    if gene_sequence[node] == sequence:
                        neighbour_list = []
                        for neighbor in G.neighbors(node):
                            if gene_sequence[neighbor] == sequence:
                                neighbour_list.append(neighbor)
                        if len(neighbour_list) == 1:
                            end_list.append(node)
            print('Constructed list of ends for this gene')
            return end_list
        except KeyError:
            raise KeyError('Given gene is not present')
            
    def order_sequences(self,start_gene,graph, gene_present_dict, gene_sequence):
        end_list = self.ends_of_sequence(start_gene,graph, gene_present_dict, gene_sequence)
        
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
                    if child in gene_present_dict:
                        if gene_present_dict[child] and gene_sequence[child]!=gene_sequence[gene]:
                            yield child
                            break
                        else:
                            yield parent, child
                            visited.append(child)
                            queue.append((child, neighbors(child)))
            except StopIteration:
                queue.popleft()
        
    
    def closest_gene(self,start_gene, graph, gene_present_dict, gene_sequence):
        output_list = self.order_sequences(start_gene, graph, gene_present_dict, gene_sequence)
        
        x = output_list.__next__()
        while type(x) == tuple:
            x = output_list.__next__()
        print('Found closest gene')
        return x
    
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
        for gene in gene_sequence:
            sequence = gene_sequence[gene]
            if sequence not in sequence_list:
                sequence_list.append(sequence)
        print('Created list of all sequences present')
    
    #Construct list of unused sequences
        unused_sequences = sequence_list
        print('Constructed list of unused sequences')
    
    #Construct subgraph of all the genes that are present
        h = G.subgraph([gene for gene in gene_present_dict if gene_present_dict[gene]])
        print('Constructed subgraph')
    
    #Construct the dictionary of pairs of closest genes
        closest_genes_dict = {}
        for gene in gene_present_dict:
            if gene_present_dict[gene]:
                closest_genes_dict[end_list[0]] = self.closest_gene(end_list[0],h, gene_present_dict, gene_sequence)
                closest_genes_dict[end_list[1]] = self.closest_gene(end_list[1],h, gene_present_dict, gene_sequence)
        
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
                        other_end = [i for i in self.ends_of_sequence(end_gene, h, gene_present_dict, gene_sequence) if i!= end_gene][0]
                        sequence_2 = gene_sequence[other_end]
                        vistied_genes.append(end_gene)
                        visited_genes.append(other_end)
                        if sequence_1 == sequence_2:
                            if nx.shortest_path_length(G, end_gene, closest_genes_dict[end_gene]) > nx.shortest_path_length(G, other_end, closest_genes_dict[other_end]):
                                closest_genes_dict[end_gene] = None
                            else:
                                closest_genes_dict[other_end] = None
        print('Constructed closest_genes_dict')
    
    #Construct list of genes present in the file- only important if there are zero or one entries
        sequences_genes_present = []
        for gene in gene_present_dict:
            if gene_present_dict[gene]:
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
            k = G.subgraph([gene for gene in gene_present_dict if gene_sequence[gene]==sequence_present])
            nx.drawing.nx_pydot.write_dot(k,self.outputgraphfile)
            print('One sequence with genes present')
        else:
            if len(closest_genes_dict) == 0: #Understand why it is possible to have no elements in closest_genes_dict
                if len(G.nodes()) != 0:
                    example_present_gene = G.nodes()[0]
                    sequence = gene_sequence[example_present_gene]
                    k = G.subgraph([gene for gene in index_file if gene_sequence[gene]==sequence])
            else:
                for key in list(closest_genes_dict.keys()):
                    if closest_genes_dict[key] == None:
                        start_end = key
                        break
                print('Found starting point')
                
                unused_sequences = sequence_list
                for gene in gene_present_dict:
                    if gene_present[gene]:
                        if gene_sequence[gene] in unused_sequences:
                            unused_sequences.remove(gene_sequence[gene])
                print('Found unused sequences')
                
                for key in list(closest_genes_dict.keys()):
                    if closest_genes_dict[key]!= None:
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
                        