from create_subproblems import *
from path import *

G = toy_network_2()
tm = generate_uniform_tm(G)
num_clusters = int(np.sqrt(len(G.nodes)))
# G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict, hash_for_clusterid, meta_to_virt_dict, virt_to_meta_dict, _ = construct_subproblems(G, tm, num_clusters=num_clusters)

vis_graph(G, tm, num_clusters)

