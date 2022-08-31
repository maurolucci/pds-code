from powersimdata import Grid
import networkx as nx
import numpy as np
import sys

if __name__ == '__main__':
    interconnects = ['Eastern', 'Western', 'Texas', 'USA']
    if len(sys.argv) > 1:
        interconnects = sys.argv[1:]
    for con in interconnects:
        grid = Grid(con)
        graph = nx.Graph()
        graph.add_edges_from(np.array(grid.branch[['from_bus_id', 'to_bus_id']]))
        for n in graph.nodes:
            graph.nodes[n]['zero_injection'] = 1
        for plant in grid.plant['bus_id']:
            graph.nodes[plant]['zero_injection'] = 0
        for n in graph.nodes():
            if grid.bus['Pd'][n] > 0 or grid.bus['Qd'][n] > 0:
                graph.nodes[n]['zero_injection'] = 0
        for line in nx.generate_graphml(graph):
            print(line)
        #print(len(grid.bus), len(graph))
        #print(len(grid.branch), graph.num_edges())
        #print(len[n for n in graph.nodes() if graph.nodes[n]['zero_injection']])
