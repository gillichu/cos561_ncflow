import glob
from collections import defaultdict 
import numpy as np
from ncflow_singleiter import solve
# from ncflow_singleiter_alt import solve
from preprocess import * 
from benchmarks import preprocess_tz_files
import time 

def solver_wrapper(G, tm, max_num_iters, num_clusters):
    current_total_demand = 0
    current_commodity_demands = defaultdict(dict)
    orig_G  = G.copy()
    current_G = G
    current_tm = tm
    num_nodes = len(G.nodes())
    #num_clusters = int(np.sqrt(num_nodes))
    # num_clusters = 60

    total_demand_fulfilled = 0 
    start_time = time.time()

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
        seen = {}
        for commod in sdict.keys():
            path_alloc = sdict[commod]
            print(commod, path_alloc, [current_G[u1][u2]['capacity'] for ((u1, u2), path_flow) in path_alloc])
            if len(path_alloc) > 0 : 
                # demand_fulfilled = min([min(current_tm[u1][u2],path_flow,current_G[u1][u2]['capacity']) for ((u1, u2), path_flow) in path_alloc])
                for edge in path_alloc:
                    #print(edge)]
                    ((u1, u2), path_flow) = edge
                    if current_G[u1][u2]['capacity'] >= path_flow:
                        current_G[u1][u2]['capacity'] -= path_flow
                    else:
                        current_G[u1][u2]['capacity'] = 0
                    if (u1, u2, orig_G[u1][u2]['capacity']) not in seen:
                        seen[(u1,u2, orig_G[u1][u2]['capacity'])] = 0.
                    seen[(u1,u2, orig_G[u1][u2]['capacity'])] += path_flow
                    # demand_fulfilled = path_alloc[-1][-1]
                    # if demand_fulfilled > current_tm[u1][u2]:
                    #     print(current_tm[u1][u2], demand_fulfilled)
                    #     raise Exception('problem')
                # demand_fulfilled = min([min(current_tm[u1][u2],path_flow) for ((u1, u2), path_flow) in path_alloc])
                
                # update the current commodity demands
                #print(path_alloc)
                sources, destinations = [edge[0][0] for edge in path_alloc], [edge[0][1] for edge in path_alloc]
                for i in range(len(sources)):
                    u1 = sources[i]
                    for j in range(i, len(destinations)):
                        u2 = destinations[j]
                    demand_fulfilled = path_alloc[-1][-1]
                    
                    if current_tm[u1][u2] >= demand_fulfilled:
                        current_tm[u1][u2] -= demand_fulfilled
                        total_demand_fulfilled += demand_fulfilled
                    else:
                        total_demand_fulfilled += current_tm[u1][u2]
                        current_tm[u1][u2] = 0

        # --------------------------------------------------------------------------------
        # 3. check to see if we can quit early
        # --------------------------------------------------------------------------------
        curr_total_capacity = 0.0
        curr_total_capacity = sum(cap for _, _, cap in current_G.edges.data('capacity'))
        
        print('found total demand = ', total_demand_fulfilled)

        if curr_total_capacity == 0.0:
            print('total residual capacity equals 0.0')
            print('found total demand = ', total_demand_fulfilled)
            print('iteration number', iteridx)
            break

    print('found total demand = ', total_demand_fulfilled)
    total_original_demand = sum(sum(tm))
    print('total orginal demand', total_original_demand)

    return total_demand_fulfilled, time.time()-start_time


if __name__ == '__main__':
    max_num_iters = 1#10

    #graphname = "cogentco"
    #graphfile = "/Users/gc3045/cos561/cos561_ncflow/topologies/Cogentco.graphml"
    #tmfile = "/Users/gc3045/cos561/ncflow/traffic-matrices/uniform/Cogentco.graphml_uniform_1022466024_16.0_0.06_traffic-matrix.pkl"

    #graphname = "cogentco"
    #graphfile = "../topologies/Cogentco.graphml"
    #tmfile = "../traffic-matrices/uniform/Cogentco.graphml_uniform_1022466024_16.0_0.06_traffic-matrix.pkl"

    graphname = "Uninett2010"
    # graphname = "Kdl"
    # graphname = "Cogentco"
    graphfile = f'/Users/simran/cos561/proj/cos561_ncflow/topologies/{graphname}.graphml'
    #tmfile = "../../ncflow/traffic-matrices/uniform/Uninett2010.graphml_uniform_1089401497_4.0_0.46_traffic-matrix.pkl"

    ### GILLIAN WILL DO BIMODAL

    dist_types = ["uniform"]#["bimodal", "gravity", "uniform"]
    for dist_type in dist_types:
        print(f'dist_type: {dist_type}')
        tmfile = glob.glob(f'/Users/simran/cos561/proj/other/traffic-matrices/{dist_type}/{graphname}*')[0]

        print(f'tmfile: {tmfile}')

        outfile = f'{graphname}_{dist_type}_results.txt'

        tm = read_traffic_matrix(tmfile)
        G = read_graphml(graphfile)

        num_nodes = len(G.nodes())
        print(f"{graphname} contains {num_nodes} nodes")

        G, tm, num_nodes = preprocess_tz_files(graphname, G, tm, num_nodes)

        num_clusters_unit = num_nodes / (9.0*2)
        num_clusters = [int(i * num_clusters_unit) for i in range(1, 6)]
        print("num_clusters", num_clusters)

        #total_demand, time_taken = solver_wrapper(G, tm, max_num_iters, 10)
        #print("total demand", total_demand)
        #print("time taken", time_taken)

        demand_results = []
        time_results = []

        for num_cluster in num_clusters:

            tm = read_traffic_matrix(tmfile)
            G = read_graphml(graphfile)
            num_nodes = len(G.nodes())
            print(f"{graphname} contains {num_nodes} nodes")

            G, tm, num_nodes = preprocess_tz_files(graphname, G, tm, num_nodes)

            print('num_cluster', num_cluster)
            total_demand, time_taken = solver_wrapper(G, tm, max_num_iters, int(num_cluster))
            demand_results.append(total_demand)
            time_results.append(time_taken)

        with open(outfile, 'w+') as f:
            f.write(str(num_clusters) + "\n")
            f.write(str(demand_results) + "\n")
            f.write(str(time_results) + "\n")

