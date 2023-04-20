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

# Return the selected edge between pairs of meta nodes with its capacity for iter_id-th iteration
# result[(u_meta,v_meta)] = (selected edge (u,v),capacity of (u,v))
def select_inter_edge(G, agg_edge_dict, iter_id):
    cap_list = nx.get_edge_attributes(G,'capacity')
    selected_inter_edge = defaultdict(dict)
    for (u_meta,v_meta) in agg_edge_dict.keys():
        edges = agg_edge_dict[(u_meta,v_meta)] 
        i = iter_id % len(edges)
        selected_inter_edge[(u_meta,v_meta)] = (edges[i],cap_list[edges[i]])
    return selected_inter_edge
    
# Check if input path is a valid path, return the min edge cap, i.e. path bottleneck capacity
def compute_agg_path_cap(G_agg,path, selected_inter_edge):
    if nx.is_path(G_agg, path) == False:
        print("not a path in compute_agg_path_cap!")
        return 0

    min = 999999
    for (u_cluster,v_cluster) in path_to_edgelist(path):
        (_,cap) = selected_inter_edge[(u_cluster,v_cluster)]
        if cap < min: min = cap
    return min

# Find selected path from u_cluster to v_cluster in the iter_id-th iteration 
def path_meta_pair(G, G_agg, u_cluster, v_cluster, agg_edge_dict, iter_id,  k = None):
    paths = path_simple(G_agg, u_cluster,v_cluster,k)
    if len(paths) < 1:
        print("no path found in path_meta_pair!")
        return ([],0)
    selected_inter_edge = select_inter_edge(G, agg_edge_dict, iter_id)
    paths = sorted(paths, key = lambda path: compute_agg_path_cap(G_agg,path,selected_inter_edge), reverse=True)
    return (paths[0],compute_agg_path_cap(G_agg,paths[0],selected_inter_edge))

# Find selected path for every pair of clusters in the iter_id-th iteration
# OUTPUT: paths[(u_cluster,v_cluster)] = 
#           (the selected path represented by a list of cluster node ids, bottleneck cap of the selected path)
def path_meta(G, G_agg, num_clusters, agg_edge_dict, iter_id, k = None):
    meta_pairs = list(combinations(range(0,num_clusters), 2))
    paths = defaultdict(dict)
    for (u_cluster, v_cluster) in meta_pairs:
        if u_cluster == v_cluster: continue
        paths[(u_cluster,v_cluster)] = path_meta_pair(G, G_agg, u_cluster, v_cluster, agg_edge_dict, iter_id, k)
        paths[(v_cluster,u_cluster)] = path_meta_pair(G, G_agg, v_cluster, u_cluster, agg_edge_dict, iter_id, k)

    return paths


def toy_network_2():
    G = nx.DiGraph()
    G.add_node(0, pos=(-3, 1))
    G.add_node(1, pos=(-2, 0))
    G.add_node(2, pos=(-3, -1))
    G.add_node(3, pos=(-1, 1))
    G.add_node(4, pos=(1, 0))
    G.add_node(5, pos=(-1, -1))
    G.add_node(6, pos=(2, 1))
    G.add_node(7, pos=(3, 0))
    G.add_node(8, pos=(2, -1))

    cap1, cap2 = 1, 2
    add_bi_edge(G, 0, 1, capacity=1)
    add_bi_edge(G, 0, 2, capacity=2)
    add_bi_edge(G, 1, 2, capacity=3)
    
    add_bi_edge(G, 3, 4, capacity=1)  # 4
    add_bi_edge(G, 3, 5, capacity=2)  # 5
    add_bi_edge(G, 4, 5, capacity=3)  # ...

    add_bi_edge(G, 6, 7, capacity=1)
    add_bi_edge(G, 6, 8, capacity=2)
    add_bi_edge(G, 7, 8, capacity=3)

    add_bi_edge(G, 0, 3, capacity=1)                                                                     
    add_bi_edge(G, 2, 5, capacity=1)
    add_bi_edge(G, 4, 7, capacity=1)  # 12

    return G


if __name__ == '__main__':
    G = toy_network_2()
    tm = generate_uniform_tm(G)
    # vis_graph(G)

    num_clusters = int(np.sqrt(len(G.nodes)))
    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict, clusters_commodities_dict = construct_subproblems(G, tm, num_clusters=num_clusters)
    vis_graph(G,orig_to_agg_node=orig_to_agg_node)
    vis_graph(G_agg)

    print('agg_to_orig_nodes ', agg_to_orig_nodes, '\n') 
    print('orig_to_agg_node ', orig_to_agg_node, '\n') 
    print('agg_commodities_dict ', agg_commodities_dict, '\n') 
    print('agg_edge_dict ', agg_edge_dict, '\n')

    selected_inter_edge = select_inter_edge(G, agg_edge_dict, 0)
    print('selected_inter_edge ', selected_inter_edge, '\n')

    paths = path_meta(G, G_agg, num_clusters, agg_edge_dict, 0)
    print('path_meta ', paths, '\n')

    print(compute_agg_path_cap(G_agg,[0,1,2],selected_inter_edge))
    print(compute_agg_path_cap(G_agg,[2,0,1],selected_inter_edge))





# # Find all paths from u_cluster to v_cluster
# def path_meta_pair(G, agg_edge_dict, agg_to_orig_nodes, u_cluster, v_cluster, k = None):
#     # First, add direct edges between u and v as length-1 paths
#     paths = [list(e) for e in agg_edge_dict[(u_cluster, v_cluster)]]

#     # Create G' by contracting nodes in u_cluster into a single node, and same with v_cluster
#     newG = G
#     u_nodes = agg_to_orig_nodes[u_cluster]
#     for i in range(1,len(u_nodes)):
#         newG = nx.contracted_nodes(newG,u_nodes[0],u_nodes[i],self_loops=False)
    
#     v_nodes = agg_to_orig_nodes[v_cluster]
#     for i in range(1,len(v_nodes)):
#         newG = nx.contracted_nodes(newG,v_nodes[0],v_nodes[i],self_loops=False)

#     # Find paths connecting u and v with length >= 2, these paths' endpoints are cluster_ids
#     # instread of real node_ids, so we will replace them with node_ids 
#     # **how to select which edge? a lot of details needed
#     long_paths = [path for path in path_simple(newG, u_nodes[0],v_nodes[0],k) if len(path) > 2]
#     for path in long_paths:
#         s = list(set(u_nodes).intersection(list(G.predecessors(path[1]))))[0]
#         t = list(set(v_nodes).intersection(list(G.successors(path[-2]))))[0]
#         path[0] = s
#         path[-1] = t
    
#     return paths+long_paths

# # Find paths for every pair of clusters
# # OUTPUT: paths[(u_cluster,v_cluster)] = list of all paths from u_cluster to v_cluster
# def path_meta(G, num_clusters, agg_edge_dict, agg_to_orig_nodes, k = None):
#     meta_pairs = list(combinations(range(0,num_clusters), 2))
#     paths = defaultdict(dict)
#     for (u_cluster, v_cluster) in meta_pairs:
#         paths[(u_cluster,v_cluster)] = path_meta_pair(G, agg_edge_dict, agg_to_orig_nodes, u_cluster, v_cluster, k)
#         paths[(v_cluster,u_cluster)] = path_meta_pair(G, agg_edge_dict, agg_to_orig_nodes, v_cluster, u_cluster, k)

#     return paths



# # sort agg_edge_dict: agg_edge_dict[(u_meta,v_meta)] is a list of edges between u_meta -> v_meta by decreasing edge cap
# def agg_edge_dict_sorted(G, agg_edge_dict):
#     cap_list = nx.get_edge_attributes(G,'capacity')
#     agg_edge_dict_sorted = defaultdict(dict)
#     for (u_meta,v_meta) in agg_edge_dict.keys():
#         edges = agg_edge_dict[(u_meta,v_meta)]
#         agg_edge_dict_sorted[(u_meta,v_meta)] = sorted(edges, key = lambda e: cap_list[e],reverse=True)
#     return agg_edge_dict_sorted