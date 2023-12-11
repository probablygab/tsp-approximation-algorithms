# cython: boundscheck=False, initializedcheck=False, nonecheck=False, cdivision=True
import cython
import networkx as nx
import math
import queue

def createGraphFromFile(path: str) -> nx.Graph:
    lines = open(path, 'r').readlines()

    # Init graph 
    G = nx.Graph()

    # Create nodes
    for line in lines[6:]:
        if line.strip() == 'EOF':
            break

        data = line.split()
        node, x, y = int(data[0]), float(data[1]), float(data[2])

        G.add_node(node, pos=(x, y))

    # Create edges (complete graph)
    edges = [
        (node1, node2, round(math.sqrt((data1['pos'][0] - data2['pos'][0])**2 + 
                                       (data1['pos'][1] - data2['pos'][1])**2)))
        for node1, data1 in G.nodes.items()
        for node2, data2 in G.nodes.items()
        if node1 != node2
    ]

    G.add_weighted_edges_from(edges)

    return G

def getLeaf(mst: nx.Graph) -> int:
    leaf = 0

    for node, degree in mst.degree:
        if degree == 1:
            leaf = node
            break

    return leaf
    
def getDistanceFromPath(graph: nx.Graph, path: list) -> int:
    totalDistance = 0

    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        totalDistance += graph[u][v]['weight']

    return totalDistance
    
def toHamiltonianPath(path: list) -> list:
    hpath   = []
    visited = set()

    for node in path:
        if node not in visited:
            hpath.append(node)
            visited.add(node)

    # Add return edge
    hpath.append(path[0])

    return hpath
    
def dfsIter(graph: nx.Graph, start: int):
    path    = []
    visited = set()
    stack   = [(start, graph.neighbors(start))]
    
    while stack:
        node, neighbors = stack[-1]

        # Add node to path. Search started from this node
        # Or is going back to this node (one neighbor was popped)
        path.append(node)
        visited.add(node)

        for nbor in neighbors:
            if nbor not in visited:
                stack.append((nbor, graph.neighbors(nbor)))
                break
        else:
            stack.pop()

    return path

def twiceAroundTheTree(graph: nx.Graph) -> int:
    mst = nx.minimum_spanning_tree(graph)

    # Find a leaf node (any)
    root = getLeaf(mst)

    # Get a path that goes twice around the tree
    path = dfsIter(mst, root)

    # Make shortcuts and get a hamiltonian path
    path = toHamiltonianPath(path)

    # Get distance traveled
    return getDistanceFromPath(graph, path)

def christofides(graph: nx.Graph) -> int:
    mst = nx.minimum_spanning_tree(graph)

    # Find a leaf node (any)
    root = getLeaf(mst)

    # Get odd degree nodes
    oddNodes = [node for node, degree in mst.degree if degree % 2 == 1]

    # Calculate induced subgraph from oddNodes in the original graph 
    # and get a set of min matching edges
    induced  = nx.induced_subgraph(graph, oddNodes)
    matching = nx.min_weight_matching(induced)
   
    # Create multigraph from mst and matching edges
    multigraph = nx.MultiGraph(mst)
    multigraph.add_edges_from(matching)

    # Get a path that goes through every edge only once
    path = [u for u, _ in nx.eulerian_circuit(multigraph, source=root)]

    # Make shortcuts and get a hamiltonian path
    path = toHamiltonianPath(path)

    # Get distance traveled
    return getDistanceFromPath(graph, path)
    
def PyInit_tp2():
    pass