"""
	This work is a derivative of the original SweetSEQer algorithm licensed to Oliver Serang (2012):
		https://pubmed.ncbi.nlm.nih.gov/23443135/
		
	Copyright 2021 Christopher Burgoyne

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pylab as P
import numpy as np
import networkx as nx


#take in graph data from file
graph_data_file = open("gnuplot_instructions.gpi")

nodes = []
edges = []
spect_title = ""
start_node = 0.0
for line in graph_data_file:
    #print line
    if ("title" in line):
        spect_title = line.split(" ")[1].strip("\n")
    elif ("start" in line):
        start_node = float(line.split(" ")[1])
    elif (" " in line):
        edge = line.split(" ")
        edges.append((float(edge[0]), float(edge[1]), {'label': edge[2].strip("\n")}))
    else:
        nodes.append(float(line))
graph_data_file.close()

if (len(nodes) > 0):
    #build graph
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)

    P.ion()
    fig = P.figure(2)
    P.clf()
    graph = graph.to_undirected()

    # add the m/z for the lowest m/z node
    pos = nx.spring_layout(graph)
    nx.draw_networkx_edges(graph, pos, node_color="white", with_labels=False)
    edges = graph.edges()
    nx.draw_networkx_edge_labels(graph, pos, edge_labels = dict( zip(edges, [ graph[e[0]][e[1]]['label'] for e in edges ]) ) )

    # plot the starting m/z in the corner and indicate the
    # starting node with a '*'
    P.text(pos[start_node][0], pos[start_node][1], '*')

    # add the text in the corner showing the minimum m/z peak in the chain
    ax1 = P.axes()
    P.text(0, ax1.get_ylim()[1]*1.1*0.95, '* ' + str(start_node))

    P.ylim(ax1.get_ylim()[0], ax1.get_ylim()[1]*1.1)

    # turn off the axes
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    P.savefig('Glycan_%s.png' % spect_title)
