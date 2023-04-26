from collections import defaultdict 
import numpy as np
from ncflow_singleiter import solve
from preprocess import * 
from benchmarks import preprocess_tz_files

def solver_wrapper(G, tm, max_num_iters):
    current_total_demand = 0
    current_commodity_demands = defaultdict(dict)
    current_G = G
    current_tm = tm
    num_nodes = len(G.nodes())
    num_clusters = int(np.sqrt(num_nodes))

    total_demand_fulfilled = 0 

    for iteridx in range(max_num_iters):

        print("ITERATION:", iteridx)
        # more or less ignoring what they wrote in their code btw....

        # from their paper: 

        # A simple way to do this would be to deduct the allocated flow
        # and use the residual capacity on edges in the next iteration.
        # Also, we pick different edges between clusters and/or different
        # paths on the aggregated graph in different iterations

        # print the current allocated demand
        
        # --------------------------------------------------------------------------------
        # 1. now run the different LPs
        # --------------------------------------------------------------------------------
        sdict = solve(current_G, current_tm, iteridx, num_clusters)

        # sdict[k] = kth commodity? -> path with each allocation of flow on each edge
        # sdict[(commodity_k, (src, sink, demand_k))] = [(u1, u2, path_flow), (u3, u4, path_flow), ...] 

        # --------------------------------------------------------------------------------
        # 2. update current_edge_allocation and current_commodity_demands
        # --------------------------------------------------------------------------------

        # subtract what the previous iteration allocated from both
        
        # update the capacities of G based on what flow was pushed in each edge
        for commod in sdict.keys():
            path_alloc = sdict[commod]
            if len(path_alloc) > 0 : 
                for edge in path_alloc:
                    #print(edge)
                    ((u1, u2), path_flow) = edge
                    if current_G[u1][u2]['cap'] >= path_flow:
                        current_G[u1][u2]['cap'] -= path_flow
                    else:
                        current_G[u1][u2]['cap'] = 0
                # update the current commodity demands
                #print(path_alloc)
                demand_fulfilled = path_alloc[-1][-1]
                if tm[u1][u2] >= demand_fulfilled:
                    tm[u1][u2] -= demand_fulfilled
                    total_demand_fulfilled += demand_fulfilled
                else:
                    tm[u1][u2] = 0

        # --------------------------------------------------------------------------------
        # 3. check to see if we can quit early
        # --------------------------------------------------------------------------------
        curr_total_capacity = 0.0
        curr_total_capacity = sum(cap for _, _, cap in G.edges.data('cap'))

        if curr_total_capacity == 0.0:
            print('total residual capacity equals 0.0')
            break



if __name__ == '__main__':
    max_num_iters = 6

    #graphname = "cogentco"
    #graphfile = "/Users/gc3045/cos561/cos561_ncflow/topologies/Cogentco.graphml"
    #tmfile = "/Users/gc3045/cos561/ncflow/traffic-matrices/uniform/Cogentco.graphml_uniform_1022466024_16.0_0.06_traffic-matrix.pkl"

    graphname = "uninett2010"
    graphfile = "/Users/gc3045/cos561/cos561_ncflow/topologies/Uninett2010.graphml"
    tmfile = "/Users/gc3045/cos561/ncflow/traffic-matrices/uniform/Uninett2010.graphml_uniform_1089401497_4.0_0.46_traffic-matrix.pkl"


    tm = read_traffic_matrix(tmfile)
    G = read_graphml(graphfile)

    num_nodes = len(G.nodes())

    G, tm, num_nodes = preprocess_tz_files(graphname, G, tm, num_nodes)

    solver_wrapper(G, tm, max_num_iters)
