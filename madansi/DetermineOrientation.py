class DetermineOrientation(object):
    def __init__(self,contig, visited,contig_orientation):
        
    
    def determine_orientation(self, visited, contig):
        
        if len(self.graph.neighbors(contig)) == 2:
            self.contig_orientation[contig] = self.determine_orientation_degree_2(visited, contig)  
        elif len(self.graph.neighbors(contig)) == 1:
            self.contig_orientation[contig] = self.determine_orientation_end_contig(contig)
        return self.contig_orientation[contig]

    def determine_orientation_degree_2(self, visited, contig):
        for neighbour in self.graph.neighbors(contig):
            if neighbour in visited:
                visited_contig = neighbour
            else:
                unvisited_contig = neighbour
        if self.contig_ends[contig][visited_contig][0] <= self.contig_ends[contig][unvisited_contig][0]:
            contig_orientation = 1
        else:
            contig_orientation = -1
        return contig_orientation
    
    def determine_orientation_start_contig(self, contig):
        neighbour = self.graph.neighbors(contig)
        if self.contig_ends[contig][neighbour[0]][0] >= self.contig_ends[contig][neighbour[0]][1]:
            contig_orientation = 1
        else:
            contig_orientation = -1
        return contig_orientation
    
    def determine_orientation_end_contig(self,contig):
        neighbour = self.graph.neighbors(contig)
        if self.contig_ends[contig][neighbour[0]][1] >= self.contig_ends[contig][neighbour[0]][0]:
            contig_orientation = 1
        else:
            contig_orientation = -1
        return contig_orientation