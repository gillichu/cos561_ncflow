import networkx as nx
from create_subproblems import *
from itertools import combinations
from collections import defaultdict
from itertools import tee, combinations, product
from traffic_matrix import *

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
def compute_agg_path_cap(G_agg,path, bundled_cap):
    if nx.is_path(G_agg, path) == False:
        print("not a path in compute_agg_path_cap!")
        return 0

    min = 999999
    for (u_cluster,v_cluster) in path_to_edgelist(path):
        cap = bundled_cap[(u_cluster,v_cluster)]
        if cap < min: min = cap
    return min

# Find selected path from u_cluster to v_cluster in the iter_id-th iteration 
def path_meta_pair(G, G_agg, u_cluster, v_cluster,bundled_cap, iter_id, k = None):
    paths = path_simple(G_agg, u_cluster,v_cluster,k)
    if len(paths) < 1:
        print("no path found in path_meta_pair!")
        return ([],0)
    # selected_inter_edge = select_inter_edge(G, agg_edge_dict, iter_id)
    paths = sorted(paths, key = lambda path: compute_agg_path_cap(G_agg,path,bundled_cap), reverse=True)
    path = paths[iter_id % len(paths)]
    return (path,compute_agg_path_cap(G_agg,path,bundled_cap))

# Find selected path for every pair of clusters in the iter_id-th iteration
# OUTPUT: paths[(u_cluster,v_cluster)] = 
#           (the selected path represented by a list of cluster node ids, bottleneck cap of the selected path)
def path_meta(G, G_agg, num_clusters, bundled_cap, iter_id, k = None):
    meta_pairs = list(combinations(range(0,num_clusters), 2))
    paths = defaultdict(dict)
    for (u_cluster, v_cluster) in meta_pairs:
        if u_cluster == v_cluster: continue
        paths[(u_cluster,v_cluster)] = path_meta_pair(G, G_agg, u_cluster, v_cluster, bundled_cap, iter_id, k)
        paths[(v_cluster,u_cluster)] = path_meta_pair(G, G_agg, v_cluster, u_cluster, bundled_cap, iter_id, k)

    return paths

# Find all paths for every pair of internal nodes in metanode(meta_id)
# OUTPUT: paths[(u,v)] = [path1,path2,...], where path1=[u,n1,n2,...,v]
def path_r2(meta_id,subgraph,all_v_hat_in, all_v_hat_out, agg_to_orig_nodes):
    v_hat_in, v_hat_out = all_v_hat_in[meta_id], all_v_hat_out[meta_id]
    nodes = agg_to_orig_nodes[meta_id]
    from_nodes = set(nodes).union(v_hat_in)
    to_nodes = set(nodes).union(v_hat_out)

    paths = defaultdict(dict)
    for s, t in product(from_nodes, to_nodes):
        try:
            path = path_simple(subgraph,s,t)
            paths[(s,t)] = path
        except:
            continue
    return paths

# FOR R2: compute the v_hat_ins and v_hat_outs dictionary
# e.g. for cluster id i, v_hat_ins[i] gives a list of cluster nodes that are predecessors to i
def v_hat_dict(G_agg, meta_to_virt_dict):
    v_hat_ins = defaultdict(list)
    v_hat_outs = defaultdict(list)
    for (u_meta, v_meta) in G_agg.edges:
        v_hat_ins[v_meta].append(meta_to_virt_dict[u_meta][0])
        v_hat_outs[u_meta].append(meta_to_virt_dict[v_meta][1])
    return v_hat_ins, v_hat_outs

EPS = 1e-5

def extract_sol_as_mat2(model, G, path_id_to_commod_id, all_paths):
        edge_idx = {edge: e for e, edge in enumerate(G.edges)}
        sol_mat = np.zeros(
            (len(edge_idx), len(set(path_id_to_commod_id.values()))),
            dtype=np.float32)
        for var in model.getVars():
            if var.varName.startswith('f[') and var.x > EPS:
                match = re.match(r'f\[(\d+)\]', var.varName)
                p = int(match.group(1))
                k = path_id_to_commod_id[p]
                for edge in path_to_edge_list(all_paths[p]):
                    sol_mat[edge_idx[edge], k] += var.x

        return sol_mat

def get_in_and_out_neighbors(flow_list, curr_meta_node):
	# input
	# flow_seq -> list of edges, flow allocation on those edges
	# current_meta_node

	# return set of in_neighbors, out_neighbors

	in_neighbors, out_neighbors = set(), set()
	for (u, v), l in flow_list:
		if u == curr_meta_node:
			out_neighbors.add(v)
		elif v == curr_meta_node:
			in_neighbors.add(u)
	return in_neighbors, out_neighbors

def path_to_edge_list(path):
	a, b = tee(path)
	next(b, None)
	return zip(a, b)

# Compute in flow (edge_idx => -1) or out flow (edge_idx => 0)
# to/from a given set of nodes (node_set) for a flow list
def compute_in_or_out_flow(flow_list, edge_idx, node_set={}):
    flow = 0.0
    for edge, l in flow_list:
        if edge[edge_idx] in node_set:
            flow += l
        elif edge[1 - edge_idx] in node_set:
            flow -= l

    return flow

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

def toy_network_3():
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

    add_bi_edge(G, 0, 1, capacity=1)
    add_bi_edge(G, 0, 2, capacity=2)
    add_bi_edge(G, 1, 2, capacity=3)
    
    add_bi_edge(G, 3, 4, capacity=4)  # 4
    add_bi_edge(G, 3, 5, capacity=5)  # 5
    add_bi_edge(G, 4, 5, capacity=6)  # ...

    add_bi_edge(G, 6, 7, capacity=7)
    add_bi_edge(G, 6, 8, capacity=8)
    add_bi_edge(G, 7, 8, capacity=9)

    add_bi_edge(G, 0, 3, capacity=1)                                                                     
    add_bi_edge(G, 2, 5, capacity=2)
    add_bi_edge(G, 4, 7, capacity=3)  # 12
    return G

def toy_network_4():
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

    add_bi_edge(G, 0, 1, capacity=4)
    add_bi_edge(G, 0, 2, capacity=5)
    add_bi_edge(G, 1, 2, capacity=6)
    
    add_bi_edge(G, 3, 4, capacity=7)  # 4
    add_bi_edge(G, 3, 5, capacity=8)  # 5
    add_bi_edge(G, 4, 5, capacity=9)  # ...

    add_bi_edge(G, 6, 7, capacity=10)
    add_bi_edge(G, 6, 8, capacity=11)
    add_bi_edge(G, 7, 8, capacity=12)

    add_bi_edge(G, 0, 3, capacity=4)                                                                     
    add_bi_edge(G, 2, 5, capacity=5)
    add_bi_edge(G, 4, 7, capacity=6)  # 12

    # add_bi_edge(G, 0, 6, capacity=5) 
    return G




if __name__ == '__main__':
    G = toy_network_3()
    tm = generate_uniform_tm(G)

    num_clusters = int(np.sqrt(len(G.nodes)))
    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict, clusters_commodities_dict, hash_for_clusterid, meta_to_virt_dict, virt_to_meta_dict = construct_subproblems(G, tm, num_clusters=num_clusters)
    vis_graph(G,orig_to_agg_node=orig_to_agg_node)
    # vis_graph(G_agg)

    # print('agg_to_orig_nodes ', agg_to_orig_nodes, '\n') 
    # print('orig_to_agg_node ', orig_to_agg_node, '\n') 
    # print('agg_commodities_dict ', agg_commodities_dict, '\n') 
    # print('agg_edge_dict ', agg_edge_dict, '\n')

    # selected_inter_edge = select_inter_edge(G, agg_edge_dict, 0)
    # print('selected_inter_edge ', selected_inter_edge, '\n')

    # bundle_cap = bundle_cap(G, agg_edge_dict)

    # paths = path_meta(G, G_agg, num_clusters, bundle_cap, 0)
    # print('path_meta ', paths, '\n')

    # r2_paths = path_r2(2,G_clusters_dict)
    # print('r2_paths for cluster 2: ', r2_paths, '\n')

    # v_hat_ins, v_hat_outs = v_hat_dict(G_agg)
    # print('v_hat_ins: ', v_hat_ins, '\n')
    # print('v_hat_outs: ', v_hat_outs, '\n')
    




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
