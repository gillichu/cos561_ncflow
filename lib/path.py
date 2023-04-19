import networkx as nx
from create_subproblems import *
from itertools import combinations
from collections import defaultdict

# Convert path to simple path by removing edges between duplicate nodes
def to_simple_path(path):
    while (len(set(path)) != len(path)):
        dup = [x for x in path if path.count(x) > 1][0]  
        indices = [i for i,d in enumerate(path) if d==dup]    
        del path[indices[0]:indices[-1]]
    return path

# Return the list of edges (u,v) in a path
def path_to_edgelist(path):
    edgelist = []
    for u, v in zip(path, path[1:]):
        edgelist.append((u, v))
    return edgelist

# Check if input path is a valid path, return the sum of edge capacity and min edge cap
def compute_path_cap(G, path):
    if nx.is_path(G, path) == False:
        print("not a path!")
        return 0
    cap_list = nx.get_edge_attributes(G,'capacity')
    cap = 0
    min = 999999
    for i in path_to_edgelist(path):
        cap += cap_list[i]
        if cap_list[i] < min: min = cap_list[i]
    return cap, min

# k shortest edge-disjoint simple paths between u and v
def path_simple(G, u, v, k = None):
    disj_paths = list(nx.edge_disjoint_paths(G, u, v))
    # Need to make sure paths are simple
    disj_paths = [to_simple_path(path) for path in disj_paths]
    disj_paths = sorted(disj_paths, key = lambda path: len(path))

    if k == None or len(disj_paths) < k:
        return disj_paths
    else:
        return disj_paths[:k]

# Find all paths from u_cluster to v_cluster
def path_meta_pair(G, agg_edge_dict, agg_to_orig_nodes, u_cluster, v_cluster, k = None):
    # First, add direct edges between u and v as length-1 paths
    paths = [list(e) for e in agg_edge_dict[(u_cluster, v_cluster)]]

    # Create G' by contracting nodes in u_cluster into a single node, and same with v_cluster
    newG = G
    u_nodes = agg_to_orig_nodes[u_cluster]
    for i in range(1,len(u_nodes)):
        newG = nx.contracted_nodes(newG,u_nodes[0],u_nodes[i],self_loops=False)
    
    v_nodes = agg_to_orig_nodes[v_cluster]
    for i in range(1,len(v_nodes)):
        newG = nx.contracted_nodes(newG,v_nodes[0],v_nodes[i],self_loops=False)

    # Find paths connecting u and v with length >= 2, these paths' endpoints are cluster_ids
    # instread of real node_ids, so we will replace them with node_ids 
    # **how to select which edge? a lot of details needed
    long_paths = [path for path in path_simple(newG, u_nodes[0],v_nodes[0],k) if len(path) > 2]
    for path in long_paths:
        s = list(set(u_nodes).intersection(list(G.predecessors(path[1]))))[0]
        t = list(set(v_nodes).intersection(list(G.successors(path[-2]))))[0]
        path[0] = s
        path[-1] = t
    
    return paths+long_paths

# Find paths for every pair of clusters
# OUTPUT: paths[(u_cluster,v_cluster)] = list of all paths from u_cluster to v_cluster
def path_meta(G, num_clusters, agg_edge_dict, agg_to_orig_nodes, k = None):
    meta_pairs = list(combinations(range(0,num_clusters), 2))
    paths = defaultdict(dict)
    for (u_cluster, v_cluster) in meta_pairs:
        paths[(u_cluster,v_cluster)] = path_meta_pair(G, agg_edge_dict, agg_to_orig_nodes, u_cluster, v_cluster, k)
        paths[(v_cluster,u_cluster)] = path_meta_pair(G, agg_edge_dict, agg_to_orig_nodes, v_cluster, u_cluster, k)

    return paths



def toy_network_2():
    G = nx.DiGraph()
    G.add_node(0, label='A0', pos=(-3, 1))
    G.add_node(1, label='A1', pos=(-2, 0))
    G.add_node(2, label='A2', pos=(-3, -1))
    G.add_node(3, label='B0', pos=(-1, 1))
    G.add_node(4, label='B1', pos=(1, 0))
    G.add_node(5, label='B2', pos=(-1, -1))
    G.add_node(6, label='C0', pos=(2, 1))
    G.add_node(7, label='C1', pos=(3, 0))
    G.add_node(8, label='C2', pos=(2, -1))

    cap1, cap2 = 1, 2
    add_bi_edge(G, 0, 1, capacity=cap1)
    add_bi_edge(G, 0, 2, capacity=cap1)
    add_bi_edge(G, 1, 2, capacity=cap1)
    
    add_bi_edge(G, 3, 4, capacity=cap1)
    add_bi_edge(G, 3, 5, capacity=cap1)
    add_bi_edge(G, 4, 5, capacity=cap1)

    add_bi_edge(G, 6, 7, capacity=cap1)
    add_bi_edge(G, 6, 8, capacity=cap1)
    add_bi_edge(G, 7, 8, capacity=cap1)

    add_bi_edge(G, 0, 3, capacity=cap1)                                                                     
    add_bi_edge(G, 2, 5, capacity=cap1)
    add_bi_edge(G, 4, 7, capacity=cap1) 

    return G


if __name__ == '__main__':
    G = toy_network_2()
    tm = generate_uniform_tm(G)
    vis_graph(G)

    num_clusters = int(np.sqrt(len(G.nodes)))
    print(num_clusters)
    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict, clusters_commodities_dict = construct_subproblems(G, tm, num_clusters=num_clusters)

    print('agg_edge_dict ', agg_edge_dict, '\n')
    print('agg_to_orig_nodes ', agg_to_orig_nodes, '\n') 
    print('orig_to_agg_node ', orig_to_agg_node, '\n') 
    print('agg_commodities_dict ', agg_commodities_dict, '\n') 

    paths = path_meta(G, num_clusters,agg_edge_dict, agg_to_orig_nodes)
    print(paths)
    # print(compute_path_cap(G,[0,1]))