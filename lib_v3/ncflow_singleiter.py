import os
import re
from LpSolver import *
from gurobipy import GRB, Model, quicksum
from path import *
from create_subproblems import *
from itertools import tee, combinations
# from r2 import *
from traffic_matrix import *

def solve(G, tm, iter_id, num_clusters): 
    def nodelist_to_edgelist(path):
        a, b = tee(path)
        next(b, None)
        return zip(a, b)


    def solve_lp():
        pass 


    def r1_lp(G, paths_dict, agg_commodities_dict, edge_to_bundlecap):
        ### expects input:
        ### paths_dict[(u_meta, v_meta)] = [([e1, e2, e3, ..., eN], mincap)]
        ### paths_dict contains meta node pairs that may not be directly connected
        print("path dictionary", paths_dict)

        r1_outfile = 'r1_out.txt'
        os.remove(r1_outfile)

        commodities = []
        commodidx_to_info = dict()
        r1_path_to_commodities = dict()
        r1_commodity_to_pathids = defaultdict(list)
        pathidx_to_edgelist = defaultdict(list)
        # holds meta edges that actually exist between clusters
        meta_edge_to_pathids = defaultdict(list)
        cap_list = nx.get_edge_attributes(G, 'capacity')
        print("cap list", cap_list)
        print("agg_commodities_dict.keys", agg_commodities_dict.keys())
        
        path_idx = 0

        for key in agg_commodities_dict.keys():
            commodity_idx, pathinfo = key
            s_k, t_k, d_k = pathinfo
            # print("looking up:", s_k, t_k)
            commodidx_to_info[commodity_idx] = key

            # get single path between each pair of meta nodes
            # should only have 1, since selected inter-cluster edges
            path_nodelist = paths_dict[(s_k, t_k)][0] 
            # original node edges
            for edge in nodelist_to_edgelist(path_nodelist):
                # print("meta-edge", edge, "demand:", d_k, "capacity", cap_list[edge])
                meta_edge_to_pathids[edge].append(path_idx)

                pathidx_to_edgelist[path_idx].append(edge)
                
            # get the path to commodity (for r3)
            r1_path_to_commodities[path_idx] = commodity_idx
            r1_commodity_to_pathids[commodity_idx].append(path_idx)
            
            commodities.append((commodity_idx, d_k, [path_idx]))
            path_idx += 1
        
        m = Model('max-flow: R1')
        # print("commodities", commodities)
        
        # create a variable for each path
        path_variables = m.addVars(path_idx, vtype=GRB.CONTINUOUS, lb=0.0, name='f')

        # set objective
        obj = quicksum(path_variables)
        m.setObjective(obj, GRB.MAXIMIZE)

        # add demand constraints
        for _, d_k, path_ids in commodities:
            # sum of all path variables for commodity k (only one) should be <= commodity k's demand (d_k)
            print("Adding commodity constraint:", _, "path ids", path_ids, " <= ", d_k)
            m.addConstr(quicksum(path_variables[p] for p in path_ids) <= d_k)

        # add meta_edge capacity constraints 
        for meta_edge in meta_edge_to_pathids.keys():

            # get all paths on this meta_edge
            path_indices = meta_edge_to_pathids[meta_edge]
            c_e = edge_to_bundlecap[meta_edge]
            print("capacity", c_e)
            # c_e = paths_dict[meta_edge][1]

            # ensure that all paths on a given meta_edge meet the meta_edge constraint
            constr_vars = [path_variables[p] for p in path_indices]
            print("Adding capacity constraints: physical meta edges", meta_edge, "uses path indices", path_indices, "<=", c_e)
            m.addConstr(quicksum(constr_vars) <= c_e)

        return LpSolver(m, None, r1_outfile), r1_path_to_commodities, pathidx_to_edgelist, commodidx_to_info, path_idx, commodities, meta_edge_to_pathids

    def get_solution_as_mat(G_meta,model, path_id_to_commod_id, paths, pathidx_to_edgelist):
        num_edges = len(paths.keys())
        num_paths = len(set(path_id_to_commod_id.values()))
        
        # set up edge to edge idx, these are original edges in the aggregated graph
        # edge_dict = dict()
        # for edge_idx, edge in enumerate(paths.keys()):
        #     edge_dict[edge] = edge_idx
        edge_dict = {
            edge: e
            for e, edge in enumerate(G_meta.edges())
        }

        sol_mat = np.zeros((num_edges, num_paths), dtype=np.float32)
        for var in model.getVars():
            # match var name back to path
            p = int(re.match(r'f\[(\d+)\]', var.varName).group(1))
            print("var", p)
            commodity_idx = path_id_to_commod_id[p]
            # from path_idx get edges
            for edge in pathidx_to_edgelist[p]:
                sol_mat[edge_dict[edge], commodity_idx] += var.x
        return sol_mat
        

    def get_solution_as_dict(model, pathidx_to_edgelist, commod_info_dict, path_id_to_commod_id):
        # outputs type dictionary:
        #   key: (commodity_id, (src, sink, demand))
        #   value: [((n1, n2), flow), ((n2, n3), flow), ...] 

        sol_dict_def = defaultdict(list)
        for var in model.getVars():
            path_idx = int(re.match(r'f\[(\d+)\]', var.varName).group(1))
            commod_id = path_id_to_commod_id[path_idx]
            
            k, (s_k, t_k, d_k) = commod_info_dict[commod_id]
            sol_dict_def[(k, (s_k, t_k, d_k))] += [(edge, var.x) for edge in pathidx_to_edgelist[path_idx]]

        sol_dict_def = dict(sol_dict_def)
        return sol_dict_def

    def r2_lp(current_meta_node, paths_dict, G, G_meta, all_v_hat_in, all_v_hat_out, intra_commods, agg_to_orig_nodes, r1_sol_dict, meta_commodity_dict, r1_sol_mat, meta_to_virt_dict, virt_to_meta_dict):
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
        
        subgraph_nodes = agg_to_orig_nodes[current_meta_node] #np.argwhere(partition_vector == current_meta_node).flatten()

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
                    if len(commod_ids)==0:
                        continue
                    total_demand = sum([d_k for _, (s_k, _, d_k) in orig_commod_list_in_k_meta if s_k == subgraph_node])

                    targets = list(all_v_hat_out)

                    meta_commod_to_multi_commod_ids[k_meta].append(len(multi_commodity_list))
                    multi_commodity_list.append(([subgraph_node], targets, total_demand, commod_ids))

            # incomer (a.k.a current metanode is a target)
            elif target_k_meta == current_meta_node:
                # loop through all subgraph nodes
                for subgraph_node in subgraph_nodes:
                    commod_ids = [k for k, (_, t_k, d_k) in orig_commod_list_in_k_meta if t_k == subgraph_node]
                    if len(commod_ids)==0:
                        continue
                    total_demand = sum([d_k for _, (_, t_k, d_k) in orig_commod_list_in_k_meta if t_k == subgraph_node])
                # total_demand = sum([d_k for _, (s_k, _, d_k) in orig_commod_list_in_k_meta if t_k == subgraph_node])

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
                # print('commod_ids',commod_ids)
                    
                multi_commodity_list.append(
                        ([meta_to_virt_dict[u][0] for u in meta_in],
                         [meta_to_virt_dict[u][1] for u in meta_out],
                         total_demand, commod_ids))

        # find the paths  
        if len(multi_commodity_list)==0:
            return None, multi_commodity_list, None, None    

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
        print(path_id_to_multi_commod_ids)


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
            print(f"Adding demand constraint {mc_id_to_path_id_to_var[multi_commod_id]} <= {demand}")
            m2.addConstr(quicksum(mc_id_to_path_id_to_var[multi_commod_id].values()) <= demand)


        # Add edge capacity constraints
        for u, v, c_e in G.edges.data('capacity'):
            if (u, v) in edge_to_path_ids:
                path_ids = edge_to_path_ids[(u, v)]
                constr_vars = [
                    var for p in path_ids
                    for var in path_id_to_commod_id_to_var[p].values()
                ]
                print(f"Adding edge capacity {constr_vars} <= {c_e}")
                m2.addConstr(quicksum(constr_vars) <= c_e)


        # Add meta-flow constraints
        meta_edge_inds = {
            edge: e
            for e, edge in enumerate(G_meta.edges())
            if edge[0] == current_meta_node or edge[-1] == current_meta_node
        }
        
        print(type(r1_sol_mat))
        for k_meta, multi_commod_ids_list in meta_commod_to_multi_commod_ids.items():
            s_k_meta, t_k_meta, _ = meta_commodity_list[k_meta][-1]
     
            if s_k_meta != current_meta_node:
                for v_hat_in in all_v_hat_in:
                    v_meta = virt_to_meta_dict[v_hat_in]
                    meta_in_flow = r1_sol_mat[meta_edge_inds[(v_meta, current_meta_node)], k_meta] 
                    if v_hat_in not in v_hat_in_paths or len(v_hat_in_paths[v_hat_in]) == 0:
                        pass
                    else:
                        constr_vars = [path_id_to_commod_id_to_var[p][multi_commod_id]
                            for p in v_hat_in_paths[v_hat_in]
                            for multi_commod_id in path_id_to_multi_commod_ids[p]
                            if multi_commod_id in multi_commod_ids_list]
                        print(f"Adding meta-in-flow constraint {constr_vars} <= {meta_in_flow}")
                        m2.addConstr(quicksum(constr_vars) <= meta_in_flow)

            if t_k_meta != current_meta_node:
                for v_hat_out in all_v_hat_out:
                    v_meta = virt_to_meta_dict[v_hat_out]
                    meta_out_flow = r1_sol_mat[meta_edge_inds[(current_meta_node, v_meta)], k_meta]
                    if v_hat_out not in v_hat_out_paths or len(v_hat_out_paths[v_hat_out]) == 0:
                        pass
                    else:
                        constr_vars = [
                            path_id_to_commod_id_to_var[p][multi_commod_id]
                            for p in v_hat_out_paths[v_hat_out]
                            for multi_commod_id in path_id_to_multi_commod_ids[p]
                            if multi_commod_id in multi_commod_ids_list]
                        print(f"Adding meta-out-flow constraint {constr_vars} <= {meta_out_flow}")
                        m2.addConstr(quicksum(constr_vars) <= meta_out_flow)
        
        return LpSolver(m2, None, r2_outfile), multi_commodity_list, path_id_to_multi_commod_ids, all_paths

    def extract_r2_solution_dict(model, agg_commodities_dict, multi_commodity_list, path_id_to_commod_id, all_paths, virt_to_meta_dict,intra_commods): 
        commod_id_to_meta_commod_id = {commod_key: meta_commod_key[0] for meta_commod_key, c_l in agg_commodities_dict.items() for commod_key, (_, _, _) in c_l}

        meta_sol_dict_def = defaultdict(list)
        intra_sol_dict_def = defaultdict(list)
        mc_id_to_path_id_to_flow = defaultdict(dict)
        active_meta_commodity_dict = {} 
        active_multi_commodity_srcs_dict = {}
        active_multi_commodity_targets_dict = {}

        r2_srcs_out_flow_lists_def = defaultdict(list)
        r2_targets_in_flow_lists_def = defaultdict(list)

        for mc_id, (srcs, targets, _, commod_list) in enumerate(multi_commodity_list):
            srcs_are_virtual = srcs[0] in virt_to_meta_dict
            targets_are_virtual = targets[0] in virt_to_meta_dict

            if targets_are_virtual and not srcs_are_virtual:
                active_multi_commodity_srcs_dict[mc_id] = tuple(commod_list)
            elif srcs_are_virtual and not targets_are_virtual:
                active_multi_commodity_targets_dict[mc_id] = tuple(commod_list)

            for k in commod_list:
                if k in commod_id_to_meta_commod_id: 
                    k_meta = commod_id_to_meta_commod_id[k]  
                    meta_commod_key = list(agg_commodities_dict.keys())[k_meta]
                    if k_meta not in active_meta_commodity_dict:
                        active_meta_commodity_dict[k_meta] = meta_commod_key
        for var in model.getVars():
            match = re.match(r'fp(\d+)_mc(\d+)', var.varName)
            p = int(match.group(1))
            mc = int(match.group(2)) 

            mc_id_to_path_id_to_flow[mc][p] = var.x

            srcs, targets, total_demand, commod_ids = multi_commodity_list[mc]
            srcs_are_virtual = srcs[0] in virt_to_meta_dict
            targets_are_virtual = targets[0] in virt_to_meta_dict 
          
            if srcs_are_virtual or targets_are_virtual:
                k_meta = commod_id_to_meta_commod_id[commod_ids[0]]
                meta_commod_key = list(agg_commodities_dict.keys())[k_meta]
                meta_sol_dict_def[meta_commod_key] += [(edge, var.x) for edge in path_to_edge_list(all_paths[p])]
            else:
                commod_key = (commod_ids[0], (srcs[0], targets[0], total_demand))
                intra_sol_dict_def[commod_key] += [(edge, var.x) for edge in path_to_edge_list(all_paths[p])]
            if srcs_are_virtual and not targets_are_virtual:
                    r2_targets_in_flow_lists_def[tuple(commod_ids)] += [(edge, var.x) for edge in path_to_edge_list(all_paths[p])]
            if targets_are_virtual and not srcs_are_virtual:
                r2_srcs_out_flow_lists_def[tuple(commod_ids)] += [(edge, var.x) for edge in path_to_edge_list(all_paths[p])]

        curr_r2_sol_dict = {commod: meta_sol_dict_def[commod] if commod in meta_sol_dict_def else [] for commod in active_meta_commodity_dict.values()}
      
        r2_srcs_out_flow_lists_dict = {commod: r2_srcs_out_flow_lists_def[commod] if commod in r2_srcs_out_flow_lists_def else [] for commod in list(active_multi_commodity_srcs_dict.values())}
        r2_targets_in_flow_lists_dict = {commod: r2_targets_in_flow_lists_def[commod] if commod in r2_targets_in_flow_lists_def else [] for commod in list(active_multi_commodity_targets_dict.values())}
        
        intra_sol_dict = {commod: intra_sol_dict_def[commod] if commod in intra_sol_dict_def else [] for commod in intra_commods}
        return curr_r2_sol_dict, r2_srcs_out_flow_lists_dict, r2_targets_in_flow_lists_dict, intra_sol_dict
# {commod: meta_sol_dict_def[commod] if commod in meta_sol_dict_def else [] for commod in active_meta_commodity_dict.values()}



    def reconciliation_lp(r2_solution_dict, u_meta, v_meta, G, agg_to_orig_nodes, orig_to_agg_node, meta_to_virt_dict):
        reconciliation_outfile = 'reconciliation_out.txt'
        try:
            os.remove(reconciliation_outfile)
        except:
            pass

        # create a subgraph to represent flows from u_meta to v_meta
        G_u_meta_v_meta = nx.DiGraph()
        nodes_in_u_meta, nodes_in_v_meta = set(), set()

        u_hat_in, v_hat_out = meta_to_virt_dict[u_meta][0], meta_to_virt_dict[v_meta][1] 
        print(u_meta, u_hat_in, v_meta, v_hat_out)

        all_u_meta_v_meta_inter_edges = [(u,v,G[u][v][CAPACITY]) for u in agg_to_orig_nodes[u_meta] \
                                         for v in G.successors(u) if orig_to_agg_node[v]==v_meta]
        interedge_to_idx = {(u,v): idx for idx, (u, v, _) in enumerate(all_u_meta_v_meta_inter_edges)}
        idx_to_interedge = {idx: (u,v) for idx, (u, v, _) in enumerate(all_u_meta_v_meta_inter_edges)}

        for (u, v, capacity) in all_u_meta_v_meta_inter_edges:
            if u not in nodes_in_u_meta:
                nodes_in_u_meta.add(u)
                G_u_meta_v_meta.add_node(u)
            if v not in nodes_in_v_meta:
                nodes_in_v_meta.add(v)
                G_u_meta_v_meta.add_node(v)
            G_u_meta_v_meta.add_edge(u, v, capacity=capacity)

        r2_u_meta_sol_dict, r2_v_meta_sol_dict = r2_solution_dict[u_meta], r2_solution_dict[v_meta]


        # only use meta commodity solution from R2 for flows from u_meta to v_meta
        shared_meta_commodities = set(r2_u_meta_sol_dict.keys()).intersection(set(r2_v_meta_sol_dict.keys()))
        
        for commod in shared_meta_commodities:
            print('commod',commod)
            print(v_hat_out, r2_u_meta_sol_dict[commod], sum([v==v_hat_out for ((u,v), _) in r2_u_meta_sol_dict[commod]]) )
            print(u_hat_in, r2_v_meta_sol_dict[commod], sum([u==u_hat_in for ((u,v), _) in r2_v_meta_sol_dict[commod]]))
            print()

        shared_meta_commodities = [commod for commod in shared_meta_commodities if sum([v==v_hat_out for ((u,v), _) in r2_u_meta_sol_dict[commod]]) > 0]
        shared_meta_commodities = [commod for commod in shared_meta_commodities if sum([u==u_hat_in for ((u,v), _) in r2_v_meta_sol_dict[commod]]) > 0]
        print(shared_meta_commodities)

        print(shared_meta_commodities)
        m = Model(f"reconciliation [meta_node (out) = {u_meta}][meta_node (in) = {v_meta}]")

        num_edges, num_commodities = len(all_u_meta_v_meta_inter_edges), len(shared_meta_commodities)
        commodity_vars = m.addVars(num_edges, num_commodities, vtype=GRB.CONTINUOUS, lb=0.0, name='f')

        # edge capacity constraints
        for edge_idx, (u, v, capacity) in enumerate(all_u_meta_v_meta_inter_edges):
            print(f"Adding capacity constraint: edge {edge_idx} <= {capacity}")
            m.addConstr(commodity_vars.sum(edge_idx, '*') <= capacity, f"capacity [edge_idx_{edge_idx}][u {u}][v {v}]")


        # flow constraints
        total_u_meta_flows, total_v_meta_flows = dict(), dict()
        for k_local, shared_meta_commodity in enumerate(shared_meta_commodities):
            u_meta_flows, v_meta_flows = r2_u_meta_sol_dict[shared_meta_commodity], r2_v_meta_sol_dict[shared_meta_commodity]        
            k_meta, (s_k, t_k, d_k) = shared_meta_commodity
            total_u_meta_flows[k_meta], total_v_meta_flows[k_meta] = 0.0, 0.0

            # sum of outflow from u <= f_{k_local}
            u_outflows_dict = defaultdict(float)
            for (u,v), outflow in u_meta_flows:
                if v==v_hat_out:
                    u_outflows_dict[u] += outflow

            for u in nodes_in_u_meta:
                outflow = max(u_outflows_dict[u],0.0)
                total_u_meta_flows[k_meta] += outflow

                outgoing_inter_edges_to_idx = [interedge_to_idx[(u,v)] for v in G_u_meta_v_meta.successors(u)]
                print(f"Adding outflow constraint: (sum of commods for edges {outgoing_inter_edges_to_idx}) <= {outflow}")
                m.addConstr(quicksum(commodity_vars[edge, k_local] for edge in outgoing_inter_edges_to_idx) <= outflow, f"outflow [u {u}]][k_meta {k_meta}]")

            # sum of inflow to v <= f_{k_local}
            v_inflows_dict = defaultdict(float)
            for (u,v), inflow in v_meta_flows:
                if u==u_hat_in: 
                    v_inflows_dict[v] += inflow

            for v in nodes_in_v_meta:
                inflow = max(v_inflows_dict[v],0.0)
                total_v_meta_flows[k_meta] += inflow

                incoming_inter_edges_to_idx = [interedge_to_idx[(u,v)] for u in G_u_meta_v_meta.predecessors(v)]
                print(f"Adding inflow constraint: (sum of commods for edges {incoming_inter_edges_to_idx}) <= {inflow}")
                m.addConstr(quicksum(commodity_vars[edge, k_local] for edge in incoming_inter_edges_to_idx) <= inflow, f"inflow [v {v}]][k_meta {k_meta}]")

        # objective is to maximize flow
        obj = quicksum(commodity_vars)
        m.setObjective(obj, GRB.MAXIMIZE)

        return LpSolver(m, None, reconciliation_outfile), G_u_meta_v_meta, shared_meta_commodities, nodes_in_u_meta, nodes_in_v_meta, all_u_meta_v_meta_inter_edges, idx_to_interedge

    def extract_reconciliation_solution_dict(G, r2_solution_dict, agg_to_orig_nodes, orig_to_agg_node, all_v_hat_in, all_v_hat_out):
        num_clusters = len(agg_to_orig_nodes)
        all_metanode_pairs = [(u_meta, v_meta) for u_meta in range(num_clusters) for v_meta in range(num_clusters) if u_meta != v_meta]

        reconciliation_solutions_dicts = {}
        for (u_meta, v_meta) in all_metanode_pairs:
            reconciliation_solver, G_u_meta_v_meta, shared_meta_commodities, nodes_in_u_meta, nodes_in_v_meta, all_u_meta_v_meta_inter_edges, idx_to_interedge = reconciliation_lp(r2_solution_dict, u_meta, v_meta, G, agg_to_orig_nodes, orig_to_agg_node,meta_to_virt_dict)
            reconciliation_solver.solve_lp(Method.BARRIER)
            
            model = reconciliation_solver._model
            l = []
            for var in model.getVars():
                match = re.match(r'{}\[(\d+),(\d+)\]'.format('f'), var.varName)
                edge_idx, k_local = int(match.group(1)), int(match.group(2))
                edge = idx_to_interedge[edge_idx] #list(G_u_meta_v_meta.edges)[edge_idx]#all_u_meta_v_meta_inter_edges[edge_idx]

                k_meta, (s_k, t_k, d_k) = shared_meta_commodities[k_local]
                l.append((edge, k_meta, s_k, t_k, d_k, var.x))
                         
            sol_dict_def = defaultdict(list)
            for edge, k_meta, s_k, t_k, d_k, flow in l:
                sol_dict_def[(k_meta, (s_k, t_k, d_k))].append((edge, flow))

            sol_dict_def = dict(sol_dict_def)
            print(u_meta, v_meta,len(shared_meta_commodities))
            reconciliation_solutions_dicts[(u_meta, v_meta)] = {commod: sol_dict_def[commod] if commod in sol_dict_def else [] for commod in shared_meta_commodities}

        # format: (u_meta, v_meta) --> (k_meta, (s_k, t_k, d_k)) --> (edge, flow)
        return reconciliation_solutions_dicts

    def kirchoffs_lp(meta_commod_key, commodity_list, r2_src_out_flows, r2_target_in_flows):
        kirchoff_outfile = 'kirchoff_out.txt'

        all_commod_ids = [commod_key[0] for commod_key in commodity_list]
        commod_id_to_ind = {k: i for i, k in enumerate(all_commod_ids)}
        _, (u_meta, v_meta, _) = meta_commod_key

        m = Model('Kirchoff\'s Law, meta-nodes {} (out) and {} (in)'.format(
            u_meta, v_meta))

        commod_vars = m.addVars(len(all_commod_ids), vtype=GRB.CONTINUOUS, lb=0.0, name='f')

        obj = quicksum(commod_vars)
        m.setObjective(obj, GRB.MAXIMIZE)

        # Demand constraints
        for k, (_, _, d_k) in commodity_list:
            m.addConstr(commod_vars[commod_id_to_ind[k]] <= d_k)

        # Kirchoff constraints
        src_out_flows = r2_src_out_flows[(u_meta, v_meta)]
        for commod_ids, total_flow in src_out_flows.items():
            total_flow = max(total_flow, 0.0)
            m.addConstr(quicksum([commod_vars[commod_id_to_ind[k]] for k in commod_ids]) <= total_flow)

        target_in_flows = r2_target_in_flows[(u_meta, v_meta)]
        for commod_ids, total_flow in target_in_flows.items():
            total_flow = max(total_flow, 0.0)
            m.addConstr(quicksum([commod_vars[commod_id_to_ind[k]] for k in commod_ids]) <= total_flow)

        return LpSolver(m, None, kirchoff_outfile)


    def extract_kirchoffs_sol(model, commodity_list):
        flow_dict = {}
        for var in model.getVars():
            match = re.match(r'f\[(\d+)\]', var.varName)
            i = int(match.group(1))
            k, _ = commodity_list[i]
            flow_dict[k] = 0.0 if var.x < EPS else var.x

        return flow_dict

    def divide_into_multi_commod_flows(multi_commod_flow_lists, src_or_target_idx, orig_to_agg_node):
        multi_commod_flows_per_k_meta = defaultdict(dict)

        for commod_ids, flow_list in multi_commod_flow_lists.items():
            k = commod_ids[0]
            rest = commodity_list[k][-1]
            src_or_target = rest[src_or_target_idx]
            src, target = rest[0], rest[1]
            u_meta, v_meta = orig_to_agg_node[src], orig_to_agg_node[target]

            multi_commod_flows_per_k_meta[(u_meta, v_meta)][commod_ids] = compute_in_or_out_flow(flow_list, src_or_target_idx, {src_or_target})
        return dict(multi_commod_flows_per_k_meta)

    def total_flow1(list):
        flow = 0
        for (_,f) in list:
            flow += f
        return flow

    def r3_lp(path_idx, commodities, meta_edge_to_pathids, reconciliation_sol_dict):

        # Setting up the model
        r3_outfile = 'r3_out.txt'
        try: os.remove(r3_outfile)
        except: pass
        m = Model('max-flow: R3')
        m.Params.LogToConsole = 0
        # print("commodities", commodities)
        
        # create a variable for each path
        path_variables = m.addVars(path_idx, vtype=GRB.CONTINUOUS, lb=0.0, name='f')

        # set objective
        obj = quicksum(path_variables)
        m.setObjective(obj, GRB.MAXIMIZE)

        # add demand constraints
        for _, d_k, path_ids in commodities:
            # sum of all path variables for commodity k (only one) should be <= commodity k's demand (d_k)
            # print("Adding commodity constraint:", _, "path ids", path_ids, " <= ", d_k)
            m.addConstr(quicksum(path_variables[p] for p in path_ids) <= d_k)

        # add meta_edge capacity constraints 
        for meta_edge in meta_edge_to_pathids.keys():

            # get all paths on this meta_edge
            path_indices = meta_edge_to_pathids[meta_edge]
            recon_flow = 0
            recon_all = reconciliation_sol_dict[meta_edge]
            for commod in recon_all:
                recon_flow += total_flow1(recon_all[commod])
            # ensure that all paths on a given meta_edge meet the meta_edge reconciliation flow
            constr_vars = [path_variables[p] for p in path_indices]
            # print("Adding capacity constraints: physical meta edges", meta_edge, "uses path indices", path_indices, "<=", c_e)
            m.addConstr(quicksum(constr_vars) <= recon_flow)

        return LpSolver(m, None, r3_outfile)

    def compute_flow_per_inter_commod(meta_commodity_dict, inter_sol_dict, kirchoff_flow_per_commod, waterfall_memoized):
        # Calculate flow per commod that traverses between two or more meta-nodes
        flow_per_inter_commod = {}
        for meta_commod_key, orig_commod_list_in_k_meta in meta_commodity_dict.items():
            commod_ids = [key[0] for key in orig_commod_list_in_k_meta]
            meta_flow_list = inter_sol_dict[meta_commod_key]
            r3_meta_flow = compute_in_or_out_flow(
                meta_flow_list, 0, {meta_commod_key[-1][0]})
            print('meta_commod_key', meta_commod_key, 'r3_meta_flow',r3_meta_flow)
            sum_of_kirchoff_flows = sum(
                kirchoff_flow_per_commod[k] for k in commod_ids)

            if abs(sum_of_kirchoff_flows - r3_meta_flow) < EPS:
                for k in commod_ids:
                    flow_per_inter_commod[k] = kirchoff_flow_per_commod[k]
            else:
                waterfall = waterfall_memoized()
                adjusted_commods = [(key[0], (key[-1][0], key[-1][1], kirchoff_flow_per_commod[key[0]]))
                                    for key in orig_commod_list_in_k_meta]

                for k, _ in adjusted_commods:
                    waterfall_flow = waterfall(r3_meta_flow, k, adjusted_commods)
                    if waterfall_flow > EPS:
                        flow_per_inter_commod[k] = waterfall_flow
                    else:
                        flow_per_inter_commod[k] = 0.0

        return flow_per_inter_commod
    def waterfall_memoized():
        # Memoize results in demand_satisfied
        demand_satisfied = {}

        def fn(flow_val, k, commods):
            if k in demand_satisfied:
                return demand_satisfied[k]

            EPS = 1e-6
            demand_remaining = {commod[0]: commod[-1][-1] for commod in commods}
            flow_remaining = flow_val
            sorted_commods = [commod[0] for commod in sorted(commods, key=lambda x: x[-1][-1])]
            while len(demand_remaining) > 0:
                k_smallest = sorted_commods[0]
                flow_to_assign = min(flow_remaining / len(commods), demand_remaining[k_smallest])
                for commod_id, (_, _, orig_demand) in commods:
                    if commod_id not in demand_remaining:
                        continue
                    demand_remaining[commod_id] -= flow_to_assign
                    if abs(demand_remaining[commod_id] - 0.0) < EPS:
                        demand_satisfied[commod_id] = orig_demand
                        del demand_remaining[commod_id]
                        sorted_commods.remove(commod_id)
                    flow_remaining -= flow_to_assign
                if abs(flow_remaining - 0.0) < EPS:
                    break
            for commod_id, (_, _, orig_demand) in commods:
                if commod_id in demand_remaining:
                    demand_satisfied[commod_id] = orig_demand - demand_remaining[commod_id]

            return demand_satisfied[k]

        return fn


#def solve(G, tm, iter_id, num_clusters): 


    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict, hash_for_clusterid, meta_to_virt_dict, virt_to_meta_dict, commodity_list = construct_subproblems(G, tm, num_clusters=num_clusters)

    #print("G clusters dict", [(k, G_clusters_dict[k].nodes()) for k in G_clusters_dict])
    print("agg_edge_dict", agg_edge_dict)
    
    # bundle capacities on inter_edges
    edge_to_bundlecap = bundle_cap(G, agg_edge_dict)
    print("edge_to_bundlecap", edge_to_bundlecap)

    # select paths for r1, this iteration
    paths = path_meta(G, G_agg, num_clusters, edge_to_bundlecap, 0)
    print("paths", paths)
    
    r1_solver, r1_path_to_commod, pathidx_to_edgelist, commodidx_to_info, r1_path_idx, r1_commodities, r1_meta_edge_to_pathids = r1_lp(G, paths, agg_commodities_dict, edge_to_bundlecap)
    print(r1_solver.solve_lp(Method.BARRIER))
    print(r1_solver._model.objVal)

    r1_sol_dict = get_solution_as_dict(r1_solver._model, pathidx_to_edgelist, commodidx_to_info, r1_path_to_commod)
    r1_sol_mat = get_solution_as_mat(G_agg, r1_solver._model, r1_path_to_commod, paths, pathidx_to_edgelist)

    inter_edges = select_inter_edge(G, agg_edge_dict, iter_id)

    all_v_hat_in, all_v_hat_out = v_hat_dict(G_agg, meta_to_virt_dict)
    for (u_meta, v_meta), ((u,v),cap) in inter_edges.items():
        G_hat_u_meta = G_clusters_dict[u_meta]
        v_hat_out = meta_to_virt_dict[v_meta][1]
        G_hat_u_meta.add_node(v_hat_out)
        G_hat_u_meta.add_edge(u, v_hat_out, capacity=cap)

        G_hat_v_meta = G_clusters_dict[v_meta]
        u_hat_in = meta_to_virt_dict[u_meta][0]
        G_hat_v_meta.add_node(u_hat_in)
        G_hat_v_meta.add_edge(u_hat_in, v, capacity=cap)
    
    
    # format: u_meta --> (k_meta, (s_k, t_k, d_k)) --> (true edge, flow)
    intra_sols_dicts = []
    r2_solution_dict = {}
    r2_srcs_out_flow_lists = defaultdict(list)
    r2_targets_in_flow_lists = defaultdict(list)
    
    r2_models = []
    r2_mc_lists = []
    r2_paths = []

    for current_meta_node in G_agg.nodes:
        paths_dict = path_r2(current_meta_node,G_clusters_dict[current_meta_node],all_v_hat_in, all_v_hat_out, agg_to_orig_nodes)
        intra_commods = clusters_commodities_dict[current_meta_node]
        v_hat_in, v_hat_out = all_v_hat_in[current_meta_node], all_v_hat_out[current_meta_node]
        r2_solver, multi_commodity_list, path_id_to_commod_id, all_paths = r2_lp(current_meta_node, paths_dict, G, G_agg, v_hat_in, v_hat_out, intra_commods, agg_to_orig_nodes, r1_sol_dict, agg_commodities_dict, r1_sol_mat, meta_to_virt_dict, virt_to_meta_dict)
        if len(multi_commodity_list)==0:
             continue 
        r2_solver.solve_lp(Method.BARRIER)
        r2_models.append(r2_solver._model)
        r2_mc_lists.append(multi_commodity_list)
        r2_paths.append(all_paths)

        model = r2_solver._model
        print('current_meta_node ',current_meta_node)
        
        curr_r2_sol_dict, r2_srcs_out_flow_lists_dict, r2_targets_in_flow_lists_dict, intra_sol_dict = extract_r2_solution_dict(model, agg_commodities_dict, multi_commodity_list, path_id_to_commod_id, all_paths, virt_to_meta_dict,intra_commods)
        r2_solution_dict[current_meta_node] =  curr_r2_sol_dict
        intra_sols_dicts.append(intra_sol_dict)

        for commod_ids, flow_list in r2_srcs_out_flow_lists_dict.items():
            r2_srcs_out_flow_lists[commod_ids] = flow_list
        for commod_ids, flow_list in r2_targets_in_flow_lists_dict.items():
            r2_targets_in_flow_lists[commod_ids] = flow_list
        print(current_meta_node, r2_solution_dict[current_meta_node])

    # format: (u_meta, v_meta) --> (k_meta, (s_k, t_k, d_k)) --> (true edge, flow)
    reconciliation_solutions_dicts = extract_reconciliation_solution_dict(G, r2_solution_dict, agg_to_orig_nodes, orig_to_agg_node,  all_v_hat_in, all_v_hat_out)
    print('Reconciliation', reconciliation_solutions_dicts)
    print('agg_to_orig_nodes',agg_to_orig_nodes)

    reconciliation_meta_out_flows = defaultdict(lambda: defaultdict(float))
    reconciliation_meta_in_flows = defaultdict(lambda: defaultdict(float))

    # meta_commod_id -> meta_node_id -> out/in flow
    for (u_meta, v_meta), reconciliation_sol_dict in reconciliation_solutions_dicts.items():
        # k_meta to flow
        after_recon_flow = dict()
        for meta_commod, flow_list in reconciliation_sol_dict.items():
            total_flow = 0.0
            for (u, v), l in flow_list:
                if orig_to_agg_node[u] == u_meta and orig_to_agg_node[v] == v_meta:
                    total_flow += l
            k_meta = meta_commod[0]
            reconciliation_meta_out_flows[k_meta][u_meta] += total_flow
            reconciliation_meta_in_flows[k_meta][v_meta] += total_flow
            after_recon_flow[k_meta] = total_flow
        print('(u_meta, v_meta)',(u_meta, v_meta))
        print('after_recon_flow',after_recon_flow)

    print('reconciliation_meta_out_flows',reconciliation_meta_out_flows)
    print('reconciliation_meta_in_flows',reconciliation_meta_in_flows)

    adjusted_meta_commodity_list = []
    r2_src_out_flows = divide_into_multi_commod_flows(r2_srcs_out_flow_lists, 0, orig_to_agg_node)
    r2_target_in_flows = divide_into_multi_commod_flows(r2_targets_in_flow_lists, 1, orig_to_agg_node)

    kirchoff_flow_per_commod = {}

    for meta_commod_key, orig_commod_list_in_k_meta in agg_commodities_dict.items():
        k_meta, (s_k_meta, t_k_meta, _) = meta_commod_key
        if len(r1_sol_dict[meta_commod_key]) == 0:
            # this seems unnecessary, but it fails an assert in _r3_lp
            adjusted_meta_commodity_list.append(
                (k_meta, (s_k_meta, t_k_meta, 0.0)))
            for k, _ in orig_commod_list_in_k_meta:
                kirchoff_flow_per_commod[k] = 0.0
            continue

        kirchoffs_solver = kirchoffs_lp(meta_commod_key, orig_commod_list_in_k_meta, r2_src_out_flows, r2_target_in_flows)
        kirchoffs_solver.solve_lp(Method.BARRIER)

        model = kirchoffs_solver._model
        flow_per_commod = extract_kirchoffs_sol(model, orig_commod_list_in_k_meta)

        adjusted_meta_demand = 0.0
        for k, _ in orig_commod_list_in_k_meta:
            adjusted_meta_demand += flow_per_commod[k]
            kirchoff_flow_per_commod[k] = flow_per_commod[k]

            adjusted_meta_commodity_list.append((k_meta, (s_k_meta, t_k_meta, adjusted_meta_demand)))

    print('adjusted_meta_commodity_list', adjusted_meta_commodity_list)
    meta_commodity_list = list(agg_commodities_dict.keys())

    r3_solver= r3_lp(r1_path_idx, r1_commodities, r1_meta_edge_to_pathids,reconciliation_solutions_dicts)
    r3_solver.solve_lp(Method.BARRIER)
    r3_sol_dict = get_solution_as_dict(r3_solver._model, pathidx_to_edgelist, commodidx_to_info, r1_path_to_commod)
    
    r3_obj_val = 0.0
    inter_sol_dict = {}
    for meta_commod_key, flow_list in r3_sol_dict.items():
        k_meta, (s_k_meta, _, _) = meta_commod_key
        flow_val = compute_in_or_out_flow(flow_list, 0, {s_k_meta})
        r3_obj_val += flow_val
        inter_sol_dict[meta_commodity_list[k_meta]] = flow_list
    

    def sol_dict():
        intra_sol_dict = {commod_key: flow_list for sol_dict in intra_sols_dicts for commod_key, flow_list in sol_dict.items()}
        commod_id_to_meta_commod_id = {commod_key: meta_commod_key[0] for meta_commod_key, c_l in agg_commodities_dict.items() for commod_key, (_, _, _) in c_l}


        sol_dict = intra_sol_dict

        flow_per_inter_commod_after_r3 = compute_flow_per_inter_commod(agg_commodities_dict,inter_sol_dict, kirchoff_flow_per_commod,waterfall_memoized)

        meta_commod_to_meta_edge_fraction_of_r3_flow = defaultdict(
            lambda: defaultdict(float))
        for meta_commod_key, orig_commod_list_in_k_meta in agg_commodities_dict.items():
            meta_flow_list = inter_sol_dict[meta_commod_key]
            r3_meta_flow = compute_in_or_out_flow(
                meta_flow_list, 0, {meta_commod_key[-1][0]})

            for (u_meta, v_meta), meta_flow_val in meta_flow_list:
                meta_commod_to_meta_edge_fraction_of_r3_flow[meta_commod_key][(
                    u_meta, v_meta)] += 0.0 if meta_flow_val < EPS or r3_meta_flow < EPS else meta_flow_val / r3_meta_flow

        for meta_node_id, model in enumerate(r2_models):
            if model is None:
                continue
            multi_commodity_list = r2_mc_lists[meta_node_id]
            r2_all_paths = r2_paths[meta_node_id]

            commod_ids_to_meta_edge_to_path_ids = defaultdict(lambda: defaultdict(list))

            commod_ids_to_path_id_to_r2_flow = defaultdict(dict)

            for var in model.getVars():
                if not var.varName.startswith('fp') or var.x <= EPS:
                    continue
                match = re.match(r'fp(\d+)_mc(\d+)', var.varName)
                r2_path_id, mc_id = int(match.group(1)), int(match.group(2))
                srcs, targets, _, commod_ids = multi_commodity_list[mc_id]
                commod_ids = tuple(commod_ids)

                if srcs[0] not in virt_to_meta_dict and targets[0] not in virt_to_meta_dict:
                    continue

                path = r2_all_paths[r2_path_id]
                start_node, end_node = path[0], path[-1]

                # Either srcs are virtual, targets are virtual, or both
                if targets[0] in virt_to_meta_dict:
                    # Leavers and transit
                    u_meta = meta_node_id
                    v_meta = virt_to_meta_dict[end_node] if end_node in virt_to_meta_dict else meta_node_id
                else:
                    # Enterers
                    u_meta = virt_to_meta_dict[start_node] if start_node in virt_to_meta_dict else meta_node_id
                    v_meta = meta_node_id

                commod_ids_to_meta_edge_to_path_ids[commod_ids][(
                    u_meta, v_meta)].append(r2_path_id)
                commod_ids_to_path_id_to_r2_flow[commod_ids][r2_path_id] = var.x

            for commod_ids, meta_edge_path_ids_dict in commod_ids_to_meta_edge_to_path_ids.items():
                meta_commod_key = meta_commodity_list[commod_id_to_meta_commod_id[commod_ids[0]]]
                path_id_to_r2_flow = commod_ids_to_path_id_to_r2_flow[commod_ids]

                for (u_meta, v_meta), path_ids in meta_edge_path_ids_dict.items():
                    mc_r2_flow = sum(
                        path_id_to_r2_flow[path_id] for path_id in path_ids)

                    for k in commod_ids:
                        commod_key = commodity_list[k]
                        if commod_key not in sol_dict:
                            sol_dict[commod_key] = []

                        # Every commodity needs to be present in the sol_dict, so we check
                        # if mc_r2_flow == 0.0 afterwards
                        if mc_r2_flow == 0.0:
                            continue
                        for path_id in path_ids:
                            # flow per path per commod = total flow per commod * (flow per path in r2 / total flow in r2) * (meta flow on meta edge in r3 / total meta flow in r3)
                            flow_per_path_per_commod = flow_per_inter_commod_after_r3[k] * (
                                path_id_to_r2_flow[path_id] / mc_r2_flow) * meta_commod_to_meta_edge_fraction_of_r3_flow[meta_commod_key][(u_meta, v_meta)]

                            if flow_per_path_per_commod == 0.0:
                                continue
                            for u, v in path_to_edge_list(r2_all_paths[path_id]):
                                if u in virt_to_meta_dict:
                                    continue
                                if v in virt_to_meta_dict:
                                    v_meta = virt_to_meta_dict[v]
                                    v = inter_edges[(meta_node_id, v_meta)][1]
                                sol_dict[commod_key].append(((u, v), flow_per_path_per_commod))

        return sol_dict
    print('\n\n\n')
    print(sol_dict())
    return sol_dict()


    
if __name__ == '__main__':
    G = toy_network_1()
    #G = toy_network_2()
    #G = toy_network_4()

    tm = generate_uniform_tm(G)
    num_clusters = int(np.sqrt(len(G.nodes)))
    iter_id = 0
    solve(G, tm, iter_id, num_clusters)



