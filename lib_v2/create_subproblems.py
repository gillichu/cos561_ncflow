import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from itertools import permutations
import pymetis
from collections import defaultdict
import pprint

CAPACITY = 'capacity'
LABEL = 'label'
POS = 'pos'

def bundle_cap(G, agg_edge_dict):
    # agg_edge_dict [(metanodeA, metanodeB): [(a1, b1), (a2, b2)], (metanodeA, metanodeC): [(a1, c1), (a1, c3)], ...]
    outdict = dict()
    for agg_edge in agg_edge_dict:
        outdict[agg_edge] = 0
        for og_graph_edge in agg_edge_dict[agg_edge]:
            outdict[agg_edge] += G.edges[og_graph_edge]['capacity']
    return outdict

def vis_graph(G, node_label='label', edge_label='capacity', orig_to_agg_node = [], title=None):
    '''
    Visualize the graph
    (NOTE: this is taken directly from their github)
    '''
    def get_node_attrs_or_default(G, attr, default_val):
        attr_dict = nx.get_node_attributes(G, attr)
        if not attr_dict:
            attr_dict = {}
        for node in G.nodes:
            if node not in attr_dict:
                if hasattr(default_val, '__call__'):
                    attr_dict[node] = default_val(node)
                else:
                    attr_dict[node] = default_val
        return attr_dict

    def uni_rand(low=-1, high=1):
        return (high - low) * np.random.rand() + low

    def random_pos(node):
        return (uni_rand(-3, 3), uni_rand(-3, 3))

    plt.figure(figsize=(14, 8))
    pos = get_node_attrs_or_default(G, 'pos', random_pos)
    colors = get_node_attrs_or_default(G, 'color', 'yellow')
    colors = [colors[node] for node in G.nodes]
    # visualize G with nodes colored by partition
    if orig_to_agg_node != []:
        allcolors = ['yellow','red','blue','green','purple','pink','white']
        for i in range(0,len(G.nodes)):
            cluster_id = orig_to_agg_node[i]
            colors[i] = allcolors[cluster_id % len(allcolors)]
    nx.draw(G, pos, node_size=1000, node_color=colors)
    node_labels = get_node_attrs_or_default(G, node_label, str)
    nx.draw_networkx_labels(G, pos, labels=node_labels)

    edge_labels = nx.get_edge_attributes(G, edge_label)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,
            label_pos=0.8)
    
    if title:
        plt.suptitle(title)

    plt.show()


def generate_commodity_list(tm):
    '''
    INPUT: a traffic matrix
    OUTPUT: commodity list
    '''
    num_rows, num_cols = tm.shape[0], tm.shape[-1]
    commodity_list = [(u, v, tm[u,v]) for u in range(num_rows) for v in range(num_cols) if u!=v]
    return list(enumerate(commodity_list))


def get_node_info(G, node, info_str=LABEL):
    '''
    INPUT: 
     - G: graph
     - node: some node in the graph
     - info_str: information to extract about node (e.g., CAPACITY, LABEL, or POS))
     OUTPUT:
     - information extracted about node as specified by info_str
    '''
    try:
       info = G.nodes.data(info_str)[node]
    except Exception as e:
        print(e) 
        print(f'Request for {info_str} is invalid')
    return info


def partition_network(G, num_clusters=3):
    '''
    INPUT: 
     - G: graph
     - num_clusters: number of clusters
    OUTPUT: partitioned network vector (where node i maps to cluster partition_vector[i])  
    '''
    num_nodes = len(G.nodes)
    adj_matrix = np.asarray(nx.adjacency_matrix(G, weight=CAPACITY).todense(), dtype=np.int32)
    adjncy = [np.argwhere(adj_matrix[node] != 0).flatten() for node in range(num_nodes)]
    
    xadj = np.cumsum([0] + [len(elem) for elem in adjncy]) #[:-1]

    eweights = [[adj_matrix[node][neighbor] for neighbor in neighbors] for node, neighbors in enumerate(adjncy)]
    eweights = np.array([elem for elems in eweights for elem in elems])
    adjncy = np.array([elem for elems in adjncy for elem in elems])
    
    _, partition_vector = pymetis.part_graph(num_clusters, adjncy=adjncy, xadj=xadj, eweights=eweights, recursive=True, options=pymetis.Options(rtype=0))
    return np.array(partition_vector)


def construct_subproblems(G, tm, num_clusters=3):
    '''
    Input: 
      - some graph G
      - number of clusters
    Output: 
       - G_agg: aggregated graph (no edges)
       - agg_edge_dict: maps (cluster_x_id, cluster_y_id) to list of all edges between cluster x and cluster y
       - agg_to_orig_nodes: maps an aggregate node (identified by a cluster_id) to a list of all nodes in the cluster
       - orig_to_agg_node: maps node to cluster_id
       - G_clusters_dict: maps cluster_id to a subgraph representing that cluster
       - agg_commodities_dict: maps 
       - clusters_commodities_dict: maps cluster to its intra-cluster commodity list
    '''
    orig_to_agg_node = partition_network(G, num_clusters=num_clusters)

    G_agg = nx.DiGraph() 
    agg_edge_dict = defaultdict(list)
    agg_to_orig_nodes = defaultdict(list) 

    G_clusters_dict = defaultdict(nx.DiGraph)
    num_clusters = len(np.unique(orig_to_agg_node))

    # Populate G_agg (without edges)
    # Construct mappings in agg_to_orig_nodes and orig_to_agg_node
    for cluster_id in range(num_clusters):
        nodes_in_cluster = np.argwhere(orig_to_agg_node == cluster_id).flatten() 
        agg_to_orig_nodes[cluster_id] = nodes_in_cluster
        for node in nodes_in_cluster:
            orig_to_agg_node[node] = cluster_id 

        agg_pos = np.mean([get_node_info(G, node, info_str=POS) for node in nodes_in_cluster], axis=0)
        G_agg.add_node(cluster_id, label=str(cluster_id), pos=agg_pos)
        
    # Populate G_clusters_dict and agg_edge_dict
    for node, node_cluster_id in enumerate(orig_to_agg_node):
        cluster = G_clusters_dict[node_cluster_id]
        cluster.add_node(node, label=get_node_info(G, node, info_str=LABEL), pos=get_node_info(G, node, info_str=POS))

        neighbor_nodes = G.successors(node)
        for neighbor_node in neighbor_nodes:
            neighbor_node_cluster_id = orig_to_agg_node[neighbor_node]
            if node_cluster_id == neighbor_node_cluster_id:
                cluster.add_node(neighbor_node, label=get_node_info(G, neighbor_node, info_str=LABEL), pos=get_node_info(G, neighbor_node, info_str=POS))
                cluster.add_edge(node, neighbor_node, capacity=G[node][neighbor_node][CAPACITY])
            else:
                agg_edge_dict[(node_cluster_id, neighbor_node_cluster_id)].append((node, neighbor_node))

    # ***Sort the edges in agg_edge_dict[(u_meta,v_meta)] by decreasing capacity
    # ***Also add edge (u_meta,v_meta) to G_agg if there is any edge between nodes in u_meta and v_meta in G
    cap_list = nx.get_edge_attributes(G,'capacity')
    for (u_meta,v_meta) in agg_edge_dict.keys():
        edges = agg_edge_dict[(u_meta,v_meta)]
        agg_edge_dict[(u_meta,v_meta)] = sorted(edges, key = lambda e: cap_list[e], reverse=True)
        if len(edges) > 0:
            G_agg.add_edge(u_meta, v_meta)

    meta_to_virt_dict = {}
    virt_to_meta_dict = {} 
    v_hat_i = max(G.nodes) + 1
    for v_meta in G_agg.nodes():
        meta_to_virt_dict[v_meta] = (v_hat_i, v_hat_i + 1)  # in and out
        virt_to_meta_dict[v_hat_i] = v_meta
        virt_to_meta_dict[v_hat_i + 1] = v_meta
        v_hat_i += 2

    # Construct inter and intra cluster commodity dicts
    commodity_list = generate_commodity_list(tm)
    agg_commodities_dict, clusters_commodities_dict = defaultdict(list), defaultdict(list)
    for k, (s_k, t_k, d_k) in commodity_list:
        s_k_cluster_id, t_k_cluster_id = orig_to_agg_node[s_k], orig_to_agg_node[t_k]
        if  s_k_cluster_id != t_k_cluster_id:
            # agg_commodities_dict[(s_k_cluster_id, t_k_cluster_id)].append(d_k)
            agg_commodities_dict[(s_k_cluster_id, t_k_cluster_id)].append((k, (s_k, t_k, d_k)))
            
        else:
            clusters_commodities_dict[s_k_cluster_id].append((k, (s_k, t_k, d_k)))
    # agg_commodities_dict = {(k, (s_k_cluster_id, t_k_cluster_id, sum(commodities))): commodities for k, ((s_k_cluster_id, t_k_cluster_id), commodities) in enumerate(agg_commodities_dict.items())}
    agg_commodities_dict = {(k, (s_k_cluster_id, t_k_cluster_id, sum(d_k for _, (_, _, d_k) in commodities))): commodities for k, ((s_k_cluster_id, t_k_cluster_id), commodities) in enumerate(agg_commodities_dict.items())}
    
    # Create a hash function/dictionary for cluster ids
    hash_for_clusterid = defaultdict(dict)
    for i in range(0,num_clusters): hash_for_clusterid[i] = i + len(G.nodes)

    return G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict, clusters_commodities_dict, hash_for_clusterid, meta_to_virt_dict, virt_to_meta_dict


### EXAMPLE 
def add_bi_edge(G, src, dest, capacity=None):
    G.add_edge(src, dest)
    G.add_edge(dest, src)
    if capacity:
        G[src][dest]['capacity'] = capacity
        G[dest][src]['capacity'] = capacity

def toy_network_1():
    G = nx.DiGraph()
    G.add_node(0, label='0', pos=(-2, 2))
    G.add_node(1, label='1', pos=(-1, 0))
    G.add_node(2, label='2', pos=(-2, -2))
    G.add_node(3, label='3', pos=(2, 2))
    G.add_node(4, label='4', pos=(1, 0))
    G.add_node(5, label='5', pos=(2, -2))

    cap1, cap2 = 10, 1
    add_bi_edge(G, 0, 3, capacity=cap1)
    add_bi_edge(G, 0, 1, capacity=cap2)                                                                     
    add_bi_edge(G, 1, 4, capacity=cap2)
    add_bi_edge(G, 1, 2, capacity=cap2)
    add_bi_edge(G, 2, 5, capacity=cap2)
    add_bi_edge(G, 3, 4, capacity=cap2)
    add_bi_edge(G, 4, 5, capacity=cap1)


    return G

def generate_uniform_tm(G, max_demand=10, seed=0):
    '''
    Generate uniform traffic matrix for network G
    '''
    np.random.seed(seed)
    num_nodes = len(G.nodes)
    tm = np.random.rand(num_nodes, num_nodes) * max_demand
    tm = tm.astype(np.float32)
    np.fill_diagonal(tm, 0.0)
    
    return tm

if __name__ == '__main__':
    
    G = toy_network_1()
    tm = generate_uniform_tm(G)
    vis_graph(G, title="original network")

    num_clusters = int(np.sqrt(len(G.nodes)))
    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict, clusters_commodities_dict,_ = construct_subproblems(G, tm, num_clusters=num_clusters)


    vis_graph(G_agg, title="aggregated network (using metanodes)")
    print('agg_edge_dict ', agg_edge_dict, '\n')
    print('agg_to_orig_nodes ', agg_to_orig_nodes, '\n') 
    print('orig_to_agg_node ', orig_to_agg_node, '\n')       
    print('agg_commodities_dict ', agg_commodities_dict, '\n') 

    for cluster_id, cluster in sorted(G_clusters_dict.items(), key=lambda x: x[0]):
        print('cluster_id ', cluster_id)
        vis_graph(cluster, title=f"cluster {cluster_id}")
        print('agg_to_orig_nodes ', agg_to_orig_nodes[cluster_id], '\n') 
        print('clusters_commodities_dict ', clusters_commodities_dict[cluster_id], '\n') 
