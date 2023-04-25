import networkx as nx
import pickle

# https://github.com/netcontract/ncflow/blob/9879bb1c36637acea2692c81c628f6c81a375d7c/lib/problem.py#L337
def read_graphml(graphml_file):
    # input: graphml file
    # output: nx.MultiDiGraph

    # read in the graph and convert to digraph
    graph = nx.read_graphml(graphml_file)
    graph = graph.to_directed()
    graph = nx.DiGraph(graph)

    scc_graph = []
    # originial paper uses largest strongly connected component as the graph
    for subgraph in nx.strongly_connected_components(graph):
        if len(graph.subgraph(subgraph)) > len(scc_graph):
            scc_graph = graph.subgraph(subgraph)

    scc_graph = nx.convert_node_labels_to_integers(scc_graph)

    # assume unit capacity for now i guess
    for i,j in scc_graph.edges():
       scc_graph[i][j]['cap'] = 1
    
    return scc_graph 

def read_traffic_matrix(pkl_file): 
    # input: traffic matrix pkl
    # output: dictionary (not convinced that this shouldn't just be a matrix)

    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)

    # convert to dict???

    return data

def generate_traffic_matrix():
    pass


# testing some stuff, ignore for now
graphml_file = "../topologies/Cogentco.graphml"
G = read_graphml(graphml_file)
print(G)

