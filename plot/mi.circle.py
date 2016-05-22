#!/usr/bin/env python
import re
import networkx as nx
import matplotlib.pyplot as plt

#plot node pairs that have mi > 0.5 and not self-pairs

def draw_graph(graph, labels=None, graph_layout='shell',
               node_size=300, node_color='red', node_alpha=0.5,
               node_text_size=10,
               edge_color='black', edge_alpha=1.0, edge_tickness=2,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    # create networkx graph
    G=nx.Graph()
    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])
    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)
    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)
    if labels is None:
        labels = range(len(graph))
    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
                                 label_pos=edge_text_pos)
    # show graph
    plt.show()
    #plt.savefig('test.png')


f = open("../MI-out.txt")
x = []
y = []
z = []
graph = []
for line in f.readlines():
    a = re.split('\t|\n',line)
    s0 = int(a[0])
    s1 = int(a[1])
    s2 = a[2]
    #x.append(s0)
    #y.append(s1)
    #z.append(s2)
    tp1 = (s0,s1)
    if float(s2) > 0.5 and s0 != s1:
        #print s2
        graph.append(tp1)
f.close()

print graph
#graph = [(0, 1), (1, 5), (1, 7), (4, 5), (4, 8), (1, 6), (3, 7), (5, 9),
#         (2, 4), (0, 4), (2, 5), (3, 6), (8, 9),(1, 2)]

# you may name your edge labels
#labels = range(65, 65+len(graph))
draw_graph(graph)

# if edge labels is not specified, numeric labels (0, 1, 2...) will be used
#draw_graph(graph)
