import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import nxviz as nv
import itertools
from networkx.algorithms import community


# Import enriched CDR3 data
data_file_path = "C:/Users/nealp/Dropbox (Partners HealthCare)/Projects/PNOIT2-1037/TCRB sequencing and HLA typing data/neals.thesis.data"
enriched_CDR3s = pd.read_csv(data_file_path + "/enriched.CDR3s.csv")

# Make the CDR3s into a list
CDR3s = enriched_CDR3s["CDR3.amino.acid.sequence"].tolist()

# Create an empty graph
G = nx.Graph()

# Add nodes to the graph
G.add_nodes_from(CDR3s)

# Get edges
edge_df = pd.read_csv("C:/Users/nealp/Dropbox (Personal)/Extension School/Thesis/R.directory/pairs.csv")

# Convert the edges to be tuples
edges = [tuple(x) for x in edge_df[['from.cdr3', 'to.cdr3']].values]
# Add the edges to the graph
G.add_edges_from(edges)

# Remove isolates: Nodes with no edges
G.remove_nodes_from(list(nx.isolates(G)))

# =============================================================================
# # Define nodes_with_m_nbrs()
# def nodes_with_m_nbrs(Graph, m):
#     """
#     Returns all nodes in graph G that have m neighbors.
#     """
#     nodes = set()
#     
#     # Iterate over all nodes in G
#     for n in Graph.nodes():
#     
#         # Check if the number of neighbors of n matches m
#         if len(list(Graph.neighbors(n))) == m:
#         
#             # Add the node n to the set
#             nodes.add(n)
#             
#     # Return the nodes with m neighbors
#     return nodes
# 
# # Compute and print all nodes in T that have 6 neighbors
# six_nbrs = nodes_with_m_nbrs(G, 6)
# print(six_nbrs)
# =============================================================================

# Calculate the degree of every node
degrees = [len(list(G.neighbors(n))) for n in G.nodes()]

# Calculate the degree centrality
deg_cent = nx.degree_centrality(G)

# Make a histogram of the distribution of degree centrality
plt.figure()
plt.hist(list(deg_cent.values()))
plt.show()

# Make a histogram of the degrees
plt.figure()
plt.hist(degrees)
plt.show()

# Plot a scatter plot looking at the number of degrees and degree centrality
# Should be perfectly linear
plt.figure()
plt.scatter(x = degrees, y = list(deg_cent.values()))
plt.show()

# =============================================================================
# # Lets find the node with the highest degree centrality
# def find_nodes_with_highest_deg_cent(G):
# 
#     # Compute the degree centrality of G: deg_cent
#     deg_cent = nx.degree_centrality(G)
#     
#     # Compute the maximum degree centrality: max_dc
#     max_dc = max(list(deg_cent.values()))
#     
#     nodes = set()
#     
#     # Iterate over the degree centrality dictionary
#     for k, v in deg_cent.items():
#     
#         # Check if the current value has the maximum degree centrality
#         if v == max_dc:
#         
#             # Add the current node to the set of nodes
#             nodes.add(k)
#             
#     return nodes
#     
# # Find the node(s) that has the highest degree centrality in T: top_dc
# top_dc = find_nodes_with_highest_deg_cent(G)
# print(top_dc)
# =============================================================================

# =============================================================================
# # Determine if a node is in a triangular relationship
# # Define is_in_triangle() 
# from itertools import combinations
# # G is a graph, n is a node: Use CDR3 string
# def is_in_triangle(G, n):
#     """
#     Checks whether a node `n` in graph `G` is in a triangle relationship or not. 
#     
#     Returns a boolean.
#     """
#     in_triangle = False
#     
#     # Iterate over all possible triangle relationship combinations
#     for n1, n2 in combinations(G.neighbors(n), 2):
#     
#         # Check if an edge exists between n1 and n2
#         if G.has_edge(n1, n2):
#             in_triangle = True
#             break
#     return in_triangle
# 
# is_in_triangle(G, "CASSLQGYGYTF")
# =============================================================================

# =============================================================================
# # Look for open triangles in our graph
# def node_in_open_triangle(G, n):
#     """
#     Checks whether pairs of neighbors of node `n` in graph `G` are in an 'open triangle' relationship with node `n`.
#     """
#     in_open_triangle = False
#     
#     # Iterate over all possible triangle relationship combinations
#     for n1, n2 in combinations(G.neighbors(n), 2):
#     
#         # Check if n1 and n2 do NOT have an edge between them
#         if not G.has_edge(n1, n2):
#         
#             in_open_triangle = True
#             
#             break
#             
#     return in_open_triangle
# 
# # Compute the number of open triangles in T
# num_open_triangles = 0
# 
# # Iterate over all the nodes in T
# for n in G.nodes():
# 
#     # Check if the current node is in an open triangle
#     if node_in_open_triangle(G, n):
#     
#         # Increment num_open_triangles
#         num_open_triangles += 1
#         
# print(num_open_triangles)
# 
# =============================================================================

# =============================================================================
# # Look at the number of cliques of a certain size
# # Function to determine how many cliques there are of at least n CDR3s
# def maximal_cliques(G, size):
#     """
#     Finds all maximal cliques in graph `G` that are of size `size`.
#     """
#     mcs = []
#     for clique in list(nx.find_cliques(G)):
#         if len(clique) >= size:
#             mcs.append(clique)
#     return mcs
# maximal_cliques(G, 10)
# =============================================================================

# Function to get nodes of interest and neighbors that can then be plotted
# G is a graph, nodes of interest is a list of nodes you want a subgraph of

# =============================================================================
# def get_nodes_and_nbrs(G, nodes_of_interest):
#      """
#      Returns a subgraph of the graph `G` with only the `nodes_of_interest` and their neighbors.
#      """
#      nodes_to_draw = []
#      
#      # Iterate over the nodes of interest
#      for n in nodes_of_interest:
#      
#          # Append the nodes of interest to nodes_to_draw
#          nodes_to_draw.append(n)
#          
#          # Iterate over all the neighbors of node n
#          for nbr in G.neighbors(n):
#          
#              # Append the neighbors of n to nodes_to_draw
#              nodes_to_draw.append(nbr)
#              
#      return G.subgraph(nodes_to_draw)
#  
# G_sub = get_nodes_and_nbrs(G, ["CAIPGQGAYGYTF"])
#  
# # Draw the subgraph to the screen
# nx.draw(G_sub, with_labels = False)
# plt.show()
# =============================================================================

# Get all the subgraphs in our graph
graphs = list(nx.connected_component_subgraphs(G))

# Limit to grpahs with at least 100 nodes
for n in graphs:
    if n.number_of_nodes() < 100:
        graphs.remove(n)

g_sub = graphs[10]
#nx.draw(g_sub)

# Try to define communities of CDR3s in each subgraph
# Create communities using modularity
communities = community.greedy_modularity_communities(g_sub)
labels = [list(x) for x in communities]

for n in g_sub.nodes():
    group = ([(i + 1) for i, labels in enumerate(labels) if n in labels])
    g_sub.node[n]["group"] = str(group)  


# get unique groups
groups = set(nx.get_node_attributes(g_sub,'group').values())
mapping = dict(zip(sorted(groups),itertools.count()))
colors = [mapping[g_sub.node[n]['group']] for n in g_sub.nodes()]

nx.draw_random(g_sub)
nx.draw(g_sub, node_color = colors)


# =============================================================================
# pos = nx.spring_layout(graphs[10])  # compute graph layout
# nx.draw_networkx_nodes(graphs[10], pos, node_size=600, node_color=colors)
# nx.draw_networkx_edges(graphs[10], pos, alpha=0.3)
# 
# =============================================================================
