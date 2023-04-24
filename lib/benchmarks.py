from create_subproblems import partition_network
from preprocess import *
from create_subproblems import * 

### RUNNING PF4

# read in graph file
graphfile = "/Users/gc3045/cos561/cos561_ncflow/topologies/Cogentco.graphml"
tmfile = "/Users/gc3045/cos561/ncflow/traffic-matrices/uniform/Cogentco.graphml_uniform_1022466024_16.0_0.06_traffic-matrix.pkl"

tm = read_traffic_matrix(tmfile)
G = read_graphml(graphfile)
num_nodes = len(G.nodes())

# vis_graph(G)

print("start with", len(G.nodes), "nodes")
# fill in the position 
node_ids = list(G.nodes())
for node_id in node_ids:
    node = G.nodes[node_id]
    if 'Latitude' not in node or 'Longitude' not in node:
        G.remove_node(node_id)
        continue

print("left with", len(G.nodes), "nodes")
node_ids = list(G.nodes())
for node_id in node_ids:
    #print("node_id", node_id)
    node = G.nodes[node_id]
    #print("node", node)

    node['pos'] = [node['Latitude'], node['Longitude']]

    #print('node[pos]', node['pos'])

for node_id in node_ids:
    node = G.nodes[node_id]
    assert 'pos' in node

print("example node:", G.nodes[0])
num_nodes = len(G.nodes)
# perform fm partitioning
# partition_vector = partition_network(G, num_clusters=42)

# construct subproblems
print("Constructing subproblems...")
G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict, clusters_commodities_dict, hash_for_clusterid = construct_subproblems(G, tm, num_clusters=num_nodes)

print("Bundling capacity...")
edge_to_bundlecap = bundle_cap(G, agg_edge_dict)
print("Finding all paths between pairs of meta nodes...")
paths = path_meta(G, G_agg, num_nodes, edge_to_bundlecap, 0, k=4)

print(paths)
