import os
import re
from LpSolver import *
from gurobipy import GRB, Model, quicksum
from path import *
from create_subproblems import *
from itertools import tee, combinations, product
from ncflow_singleiter import *

def nodelist_to_edgelist(path):
    a, b = tee(path)
    next(b, None)
    return zip(a, b)

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

def r2_lp(current_meta_node, paths_dict, G, G_meta, all_v_hat_in, all_v_hat_out, intra_commods, partition_vector, r1_sol_dict, meta_commodity_dict):
    # current_meta_node = index of current meta node
    # G = original graph
    # G_meta = meta graph
    # path_dict = dictionary, key = path id, value = list of edges
    # all_v_hat_in = list of metanodes that pass flow in to this node
    # all_v_hat_out = list of metanodes that this metanode passes flow to
    # intra_commods = list of commodities
    # partition_vector = list of length = number of nodes, tells which nodes go to which cluster

    # list of original commodities in a tuple of (source, target, demand)

    # r1 sol
    # edge index -> commodity index

    # things that are missing -> 
    # r1.sol_dict
    # r1.sol_mat

    # find intra_commods for this meta_node? <- they do this outside the function in the original code....
    r2_outfile = 'r2_out' + str(current_meta_node) + '.txt'
    print("am i actually doing anything?")
    meta_commodity_list = list(meta_commodity_dict.keys())

    # build meta_to_virt dict and virt_to_meta
    meta_to_virt_dict = {}
    virt_to_meta_dict = {}  # mapping from v_hat id to meta-node id
    v_hat_i = 0
    for v_meta in G_meta.nodes():
        meta_to_virt_dict[v_meta] = (v_hat_i, v_hat_i + 1)  # in and out
        virt_to_meta_dict[v_hat_i] = v_meta
        virt_to_meta_dict[v_hat_i + 1] = v_meta
        v_hat_i += 2

    subgraph_nodes = np.argwhere(partition_vector == current_meta_node).flatten()

    multi_commodity_list = []
    meta_commod_to_multi_commod_ids = defaultdict(list)

    for k, (source_k, target_k, demand_k) in intra_commods:
        multi_commodity_list.append(([source_k], [target_k], demand_k, [k]))


    # assume that the solution of r1 is stored as a dictionary
    # items() of dictionary are a meta_commod_key = (index?, (source, target, demand)), and flow seq
    for meta_commod_key, flow_seq in r1_sol_dict.items():
        if len(flow_seq) == 0:
            continue
        (k_meta, (source_k_meta, target_k_meta, _)) = meta_commod_key

        orig_commod_list_in_k_meta = meta_commodity_dict[meta_commod_key]


        # go through leavers / incomers / transit nodes
        # leaver (a.k.a current metanode is a source)
        if source_k_meta == current_meta_node:
            # loop through all subgraph nodes
            for subgraph_node in subgraph_nodes:
                commod_ids = [k for k, (s_k, _, d_k) in orig_commod_list_in_k_meta if s_k == subgraph_node]
            total_demand = sum([d_k for _, (s_k, _, d_k) in orig_commod_list_in_k_meta if s_k == subgraph_node])

            targets = list(all_v_hat_out)

            meta_commod_to_multi_commod_ids[k_meta].append(len(multi_commodity_list))
            multi_commodity_list.append(([subgraph_node], targets, total_demand, commod_ids))

        # incomer (a.k.a current metanode is a target)
        elif target_k_meta == current_meta_node:
            # loop through all subgraph nodes
            for subgraph_node in subgraph_nodes:
                commod_ids = [k for k, (_, t_k, d_k) in orig_commod_list_in_k_meta if t_k == subgraph_node]
            total_demand = sum([d_k for _, (s_k, _, d_k) in orig_commod_list_in_k_meta if t_k == subgraph_node])

            sources = list(all_v_hat_in)

            meta_commod_to_multi_commod_ids[k_meta].append(len(multi_commodity_list))

            multi_commodity_list.append((sources,[subgraph_node],total_demand, commod_ids))




        # if not an incomer or leaver, than the current metanode is just a transit node
        else: # current meta node is not a source nor a target
            meta_in, meta_out = get_in_and_out_neighbors(flow_seq, current_meta_node)
            if len(meta_in) == 0:
                continue

            commod_ids = [k for k, _ in orig_commod_list_in_k_meta]
            total_demand = sum([d_k for _, (_, _, d_k) in orig_commod_list_in_k_meta])
                
            meta_commod_to_multi_commod_ids[k_meta].append(len(multi_commodity_list))
                
            multi_commodity_list.append(
                    ([meta_to_virt_dict[u][0] for u in meta_in],
                     [meta_to_virt_dict[u][1] for u in meta_out],
                     total_demand, commod_ids))


    # find the paths
    all_paths = []
    next_path_id = 0
    source_target_paths = defaultdict(list)
    v_hat_in_paths = defaultdict(list)  # source -> ([path ids])
    v_hat_out_paths = defaultdict(list)  # targret -> ([path ids])
    path_id_to_multi_commod_ids = []  # path id -> [commod ids in multi_commodity_list]
    edge_to_path_ids = defaultdict(list)  # edge -> [path ids]

    for k, (s_k_list, t_k_list, _, k_list) in enumerate(multi_commodity_list):
        for s_k, t_k in product(s_k_list, t_k_list):
            if (s_k, t_k) not in source_target_paths:
                new_paths = paths_dict[(s_k, t_k)]
                for path in new_paths:
                    all_paths.append(path)
                    path_id_to_multi_commod_ids.append([])
                    source_target_paths[(s_k, t_k)].append(next_path_id)
                    if s_k in all_v_hat_in:
                        v_hat_in_paths[s_k].append(next_path_id)
                    if t_k in all_v_hat_out:
                        v_hat_out_paths[t_k].append(next_path_id)
                    for edge in path_to_edge_list(path):
                        edge_to_path_ids[edge].append(next_path_id)

                    next_path_id += 1

            for path_id in source_target_paths[(s_k, t_k)]:
                path = all_paths[path_id]
                path_id_to_multi_commod_ids[path_id].append(k)


    edge_to_path_ids = dict(edge_to_path_ids)
    v_hat_in_paths = dict(v_hat_in_paths)
    v_hat_out_paths = dict(v_hat_out_paths)


    m2 = Model('max-flow: R2, metanode {}'.format(current_meta_node))

    # turn path id and multicommod id -> variable id
    path_id_to_commod_id_to_var = defaultdict(dict)

    # turn multicommod id + path id -> variable id
    mc_id_to_path_id_to_var = defaultdict(dict)

    # turn multicommodity -> path id
    mc_id_to_path_ids = defaultdict(list)

    all_vars = []

    for path_id, multi_commod_ids in enumerate(path_id_to_multi_commod_ids):
        for multi_commod_id in multi_commod_ids:
            gb_var = m2.addVar(vtype=GRB.CONTINUOUS,lb=0.0,name='fp{}_mc{}'.format(path_id, multi_commod_id))

            path_id_to_commod_id_to_var[path_id][multi_commod_id] = gb_var
            mc_id_to_path_ids[multi_commod_id].append(path_id)
            mc_id_to_path_id_to_var[multi_commod_id][path_id] = gb_var
            all_vars.append(gb_var)
        m2.update()

    # set the objective
    obj = quicksum(all_vars) # original paper mentions several optimization functions. This is for the "simple" base case
    m2.setObjective(obj, GRB.MAXIMIZE)


    # set the demand constraints
    for multi_commod_id, (_, _, demand, _) in enumerate(multi_commodity_list):
        m2.addConstr(quicksum(mc_id_to_path_id_to_var[multi_commod_id].values()) <= demand)


    # Add edge capacity constraints
    for u, v, c_e in G.edges.data('capacity'):
        if (u, v) in edge_to_path_ids:
            path_ids = edge_to_path_ids[(u, v)]
            constr_vars = [
                var for p in path_ids
                for var in path_id_to_commod_id_to_var[p].values()
            ]
            m2.addConstr(quicksum(constr_vars) <= c_e)


    # Add meta-flow constraints
    meta_edge_inds = {
        edge: e
        for e, edge in enumerate(G_meta.edges())
        if edge[0] == current_meta_node or edge[-1] == current_meta_node
    }
    

    for k_meta, multi_commod_ids_list in meta_commod_to_multi_commod_ids.items():
        s_k_meta, t_k_meta = meta_commodity_list[k_meta]

        # 
        if s_k_meta != current_meta_node:
            for v_hat_in in all_v_hat_in:
                v_meta = virt_to_meta_dict[v_hat_in]
                # meta_in_flow = r1_sol_mat[meta_edge_inds[(
                #     v_meta, curr_meta_node)], k_meta]
                if v_hat_in not in v_hat_in_paths or len(v_hat_in_paths[v_hat_in]) == 0:
                    # if meta_in_flow > 0.001:
                    #     print('WARN: v_hat_in{} with flow{} has no grb vars'.format(
                    #         v_hat_in, meta_in_flow))
                    pass
                else:
                    constr_vars = [path_id_to_commod_id_to_var[p][multi_commod_id]
                        for p in v_hat_in_paths[v_hat_in]
                        for multi_commod_id in path_id_to_multi_commod_ids[p]
                        if multi_commod_id in multi_commod_ids_list]
                    
                    m2.addConstr(quicksum(constr_vars) <= meta_in_flow)

        if t_k_meta != current_meta_node:
            for v_hat_out in all_v_hat_out:
                v_meta = virt_to_meta_dict[v_hat_out]
                # meta_out_flow = r1_sol_mat[meta_edge_inds[(
                #     curr_meta_node, v_meta)], k_meta]
                if v_hat_out not in v_hat_out_paths or len(v_hat_out_paths[v_hat_out]) == 0:
                    # if meta_out_flow > 0.001:
                    #     print('WARN: v_hat_out{} with flow{} has now grb vars'.format(
                    #         v_hat_out, meta_out_flow))
                    pass
                else:
                    constr_vars = [
                        path_id_to_commod_id_to_var[p][multi_commod_id]
                        for p in v_hat_out_paths[v_hat_out]
                        for multi_commod_id in path_id_to_multi_commod_ids[p]
                        if multi_commod_id in multi_commod_ids_list]

                    m2.addConstr(quicksum(constr_vars) <= meta_out_flow)
    
    return LpSolver(m2, None, r2_outfile)


if __name__ == '__main__':
    G = toy_network_2()
    tm = generate_uniform_tm(G)
    num_clusters = int(np.sqrt(len(G.nodes)))
    iter_id = 0

    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict, hash_for_clusterid = construct_subproblems(G, tm, num_clusters=num_clusters)

    # print("G clusters dict", [(k, G_clusters_dict[k].nodes()) for k in G_clusters_dict])

    # select paths for r1, this iteration
    paths = path_meta(G, G_agg, num_clusters, agg_edge_dict, 0)
    
    r1_solver, r1_path_to_commod, pathidx_to_edgelist, commodidx_to_info = r1_lp(G, paths, agg_commodities_dict)
    print(r1_solver.solve_lp(Method.BARRIER))
    print(r1_solver._model.objVal)
    #print(get_solution_as_mat(r1_solver._model, r1_path_to_commod, paths, pathidx_to_edgelist))
    # print("solution as dict", get_solution_as_dict(r1_solver._model, pathidx_to_edgelist, commodidx_to_info, r1_path_to_commod))



    # test inputs
    current_meta_node = 0
    paths_dict = path_r2(0,G_clusters_dict) # dict that takes 
    G_meta = G_agg
    
    r1_sol_dict = get_solution_as_dict(r1_solver._model, pathidx_to_edgelist, commodidx_to_info, r1_path_to_commod)
    partition_vector = orig_to_agg_node

    # find all_v_hat_in and all_v_hat_out
    all_v_hat_ins, all_v_hat_outs = v_hat_dict(G_agg)
    all_v_hat_in = all_v_hat_ins[0]
    all_v_hat_out = all_v_hat_outs[0]

    intra_commods_dict = defaultdict(list)
    meta_commodity_dict = defaultdict(list)

    commodity_list = generate_commodity_list(tm)

    for k, (s_k, t_k, d_k) in commodity_list:
        s_k_meta = partition_vector[s_k]
        t_k_meta = partition_vector[t_k]
        if s_k_meta != t_k_meta:
            meta_commodity_dict[(s_k_meta, t_k_meta)].append(
                (k, (s_k, t_k, d_k)))
        else:
            intra_commods_dict[s_k_meta].append((k, (s_k, t_k, d_k)))

    intra_commods = intra_commods_dict[0]

    r2_sol = r2_lp(current_meta_node, paths_dict, G, G_meta, all_v_hat_in, all_v_hat_out, intra_commods, partition_vector, r1_sol_dict, meta_commodity_dict)
    print(r2_sol.solve_lp(Method.BARRIER))
