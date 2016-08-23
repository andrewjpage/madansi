import networkx as nx
from madansi.Assembly import Assembly
import sys

class DetermineOrientation(object):
    def __init__(self):
        pass

    def determine_orientation_degree_2(self,contig, contig_ends, visited,contig_orientation, graph, cycle, sequences):
        previous_contig = visited[len(visited) - 1]
        for neighbour in graph.neighbors(contig):
            if neighbour != previous_contig:
                next_contig = neighbour
        if contig_ends[contig][previous_contig][0] <= sequences[contig][2]/2 and contig_ends[contig][next_contig][0] > sequences[contig][2]/2:
            return 1
        elif contig_ends[contig][previous_contig][0] > sequences[contig][2]/2 and contig_ends[contig][next_contig][0] <= sequences[contig][2]/2:
            return -1
        elif contig_ends[contig][previous_contig][1] != None:
            if contig_ends[contig][previous_contig][1] > contig_ends[contig][previous_contig][0]:
                return 1
            else:
                return -1
        elif contig_ends[contig][next_contig][1] != None:
            if contig_ends[contig][next_contig][1] > contig_ends[contig][next_contig][0]:
                return -1
            else:
                return 1
        elif contig_ends[contig][previous_contig][0] < contig_ends[contig][next_contig][0]:
            return 1
        else:
            return -1      
        
    def determine_orientation_start_contig(self,contig, contig_ends, visited,contig_orientation, graph, cycle, sequences):
        neighbour = graph.neighbors(contig)
        if contig_ends[contig][neighbour[0]][0] < sequences[contig][2]/2:
            return -1
        else:
            return 1
                
    def determine_orientation_end_contig(self,contig, contig_ends, visited,contig_orientation, graph, cycle, sequences):
        neighbour = graph.neighbors(contig)
        if contig_ends[contig][neighbour[0]][0] < sequences[contig][2]/2:
            return 1
        else:
            return -1  
    
    def determine_orientation(self,contig, contig_ends, visited,contig_orientation, graph, cycle, sequences):
        if len(graph.neighbors(contig)) == 2:
            contig_orientation[contig] = self.determine_orientation_degree_2(contig, contig_ends, visited,contig_orientation, graph, cycle, sequences)  
        elif len(graph.neighbors(contig)) == 1:
            contig_orientation[contig] = self.determine_orientation_end_contig(contig, contig_ends, visited,contig_orientation, graph, cycle, sequences)
        return contig_orientation[contig]

    def determine_initial_direction_cycle(self,contig, contig_ends, visited,contig_orientation, graph, cycle, sequences):
        neighbours = graph.neighbors(contig)
        if contig_ends[contig][neighbours[0]][0] > contig_ends[contig][neighbours[1]][0]:
            next_neighbour = neighbours[0]
        else:
            next_neighbour = neighbours[1]
        contig_orientation[contig] = 1
        return next_neighbour


                    