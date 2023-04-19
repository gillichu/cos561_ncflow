import networkx as nx
from create_subproblems import *
from itertools import combinations
from collections import defaultdict

def to_simple_path(path):
    # todo
    while (len(set(path)) != len(path)):
        dup = [x for x in path if path.count(x) > 1][0]        
    return path

def compute_cap(path):
    # todo
    return path

# k shortest edge-disjoint simple paths between u and v
def path_simple(G, u, v, k = None):
    disj_paths = list(nx.edge_disjoint_paths(G, u, v))
    # Need to make sure paths are simple
    # disj_paths = [to_simple_path(path) for path in disj_paths]
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
    # vis_graph(G)
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


if __name__ == '__main__':
    G = toy_network_1()
    tm = generate_uniform_tm(G)
    vis_graph(G)

    num_clusters = int(np.sqrt(len(G.nodes)))
    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict, clusters_commodities_dict = construct_subproblems(G, tm, num_clusters=num_clusters)

    print('agg_edge_dict ', agg_edge_dict, '\n')
    print('agg_to_orig_nodes ', agg_to_orig_nodes, '\n') 
    print('orig_to_agg_node ', orig_to_agg_node, '\n') 
    print('agg_commodities_dict ', agg_commodities_dict, '\n') 

    paths = path_meta(G, num_clusters,agg_edge_dict, agg_to_orig_nodes)
    print(paths)