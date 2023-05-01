import os
import glob
from itertools import permutations 
import pickle
import networkx as nx
#from ncflow_singleiter import * 
from ncflow_singleiter import * 
from preprocess import *
from create_subproblems import * 
from path import * 
from collections import defaultdict
import time

def preprocess_tz_files(graphname, G, tm, num_nodes): 
    dropped_nodes = []
    print("start with", len(G.nodes), "nodes", "tm of shape", tm.shape)

# fill in the position 
    node_ids = list(G.nodes())
    for node_id in node_ids:
        node = G.nodes[node_id]
        if 'Latitude' not in node or 'Longitude' not in node or 'label' not in node:
            G.remove_node(node_id)
            dropped_nodes.append(node_id)
            continue

# remap all node indices in graphfile and traffic matrix
    mapping = {node_id: nidx for nidx, node_id in enumerate(G.nodes)}
    G = nx.relabel_nodes(G, mapping)
#print(G.nodes())

    node_ids = list(G.nodes())
    for node_id in node_ids:
        node = G.nodes[node_id]
        node['pos'] = [node['Latitude'], node['Longitude']]

# rename the capacity from 'cap' to 'capacity' to fit expected
    for u, v, a in G.edges(data=True):
        # print("edge", u, v, a)
        G[u][v]['capacity'] = 100 #G[u][v]['cap']
        #G[u][v]['capacity'] = G[u][v]['cap']


#print('G.nodes.data(pos)', G.nodes.data('pos'))
    for node_id in node_ids:
        assert G.nodes.data('pos')[node_id]

    print("example node:", G.nodes[0])

# perform fm partitioning
# partition_vector = partition_network(G, num_clusters=42)

# construct subproblems
    print("Constructing subproblems...")
    scc_graph = []
    for subgraph in nx.strongly_connected_components(G):
        if len(G.subgraph(subgraph)) > len(scc_graph):
            scc_graph = G.subgraph(subgraph)
    G = scc_graph

# for all nodes not in graph
    for node in node_ids:
        if node not in G.nodes:
            dropped_nodes.append(node)

    num_nodes = len(G.nodes)
    mapping = {node_id: nidx for nidx, node_id in enumerate(G.nodes)}
    G = nx.relabel_nodes(G, mapping)

    print('dropped_nodes', dropped_nodes)
    tm = np.delete(tm, dropped_nodes, axis=0)
    tm = np.delete(tm, dropped_nodes, axis=1)

    return G, tm, num_nodes


### BEGINNING OF PROCESS TOPOLOGY ZOO GRAPHS
if __name__ == '__main__':
    ### RUNNING PF4
    #tm_type = "uniform"

    # read in graph file
    ### SMALLEST (74 nodes, should be extremely quick to run.)
    #graphname = "uninett2010"
    #graphfile = "/Users/gc3045/cos561/cos561_ncflow/topologies/Uninett2010.graphml"
    #tmfile = "/Users/gc3045/cos561/ncflow/traffic-matrices/uniform/Uninett2010.graphml_uniform_1089401497_4.0_0.46_traffic-matrix.pkl"


    ### MEDIUM TOPOLOGY 

    #graphname = "cogentco"
    #graphfile = "../topologies/Cogentco.graphml"
        
    ### LARGEST TOPOLOGY (754 nodes, takes a long time to run.)
    #graphname = "kdl"
    #graphfile = "/Users/gc3045/cos561/cos561_ncflow/topologies/Kdl.graphml"
    #tmfile = "/Users/gc3045/cos561/ncflow/traffic-matrices/uniform/Kdl.graphml_uniform_1192452225_1.0_0.005_traffic-matrix.pkl"

    # vis_graph(G)

    #dist_types = ["bimodal", "gravity", "uniform"]
    #dist_types = ["bimodal", "gravity", "poisson-high-inter", "poisson-high-intra", "uniform"]
    
    # graphname = "Cogentco" #Uninett2010"
    # graphfile = f"../topologies/{graphname}.graphml"
    graphname = "Uninett2010"
    #graphname = "Cogentco"
    graphfile = f'/Users/simran/cos561/proj/cos561_ncflow/topologies/{graphname}.graphml'


    write_output = []

    #for dist_type in dist_types:
    dist_type = "uniform"
    
    outfile = f'pf4_out/{graphname}_{dist_type}_result_lib.txt'

    start_time = time.time()

    dname = f'pf4_out/{dist_type}'
    if not os.path.exists(dname):
        os.mkdir(dname)
    #tmfile = glob.glob(f'/Users/gc3045/cos561/ncflow/traffic-matrices/{dist_type}/{graphname}*')[0]
    # tmfile = '/Users/gc3045/cos561/ncflow/traffic-matrices/uniform/Cogentco.graphml_uniform_1022466024_16.0_0.06_traffic-matrix.pkl'
    #tmfile = glob.glob(f'/Users/gc3045/cos561/ncflow/traffic-matrices/{dist_type}/{graphname}*')[0]
    tmfile = glob.glob(f'/Users/simran/cos561/proj/other/traffic-matrices/{dist_type}/{graphname}*')[0]


    pkl_out = f'{dname}/{graphname}_{dist_type}_pathsformulation4.pkl'
    solndict_out = f'{dname}/{graphname}_{dist_type}_solndict.pkl'
    pf4_outfile = f'{dname}/{graphname}_{dist_type}_pf4_model_lib.txt'

    tm = read_traffic_matrix(tmfile)
    G = read_graphml(graphfile)
    num_nodes = len(G.nodes())


    G, tm, num_nodes = preprocess_tz_files(graphname, G, tm, num_nodes)

    #print("Bundling capacity...")
    edge_to_bundlecap = dict()
    for u, v, a in G.edges(data=True):
        edge_to_bundlecap[(u, v)] = G[u][v]['capacity']

    all_paths = dict()
    for u, v in permutations(list(G.nodes), 2): 
        paths = path_simple(G, u, v, k=4)
        all_paths[(u, v)] = paths

    with open(pkl_out, 'wb') as w:
        pickle.dump(all_paths, w)

    ### Need to build agg_commodities_dict
    commodities_dict = defaultdict(list)
    commodity_list = generate_commodity_list(tm)
    for k, (s_k, t_k, d_k) in commodity_list:
        commodities_dict[(k, (s_k, t_k, d_k))].append(d_k)

    ### Need to build partition_vector
    partition_vector = []
    for nidx, node in enumerate(G.nodes()):
        partition_vector.append(nidx)

    # r1_solver, r1_path_to_commod, pathidx_to_edgelist, commodidx_to_info = r1_lp(G, all_paths, commodities_dict, edge_to_bundlecap, pf4_outfile)
    r1_solver, r1_path_to_commod, pathidx_to_edgelist, commodidx_to_info, _, _, _ = r1_lp(G, all_paths, commodities_dict, edge_to_bundlecap, pf4_outfile)


    r1_solver.solve_lp(Method.BARRIER)
    objval = r1_solver._model.objVal
    
    solndict = get_solution_as_dict(r1_solver._model, pathidx_to_edgelist, commodidx_to_info, r1_path_to_commod)

    with open(solndict_out, "wb") as w:
        pickle.dump(solndict, w)

    runtime = time.time() - start_time

    write_output.append(f'{dist_type}, {objval}, {runtime}\n')

    with open(outfile, 'w+') as w:
        for s in write_output:
            w.write(s)


