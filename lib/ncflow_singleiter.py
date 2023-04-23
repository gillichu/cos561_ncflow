import os
import re
from LpSolver import *
from gurobipy import GRB, Model, quicksum
from path import *
from create_subproblems import *
from itertools import tee, combinations

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

    return LpSolver(m, None, r1_outfile), r1_path_to_commodities, pathidx_to_edgelist, commodidx_to_info

def get_solution_as_mat(model, path_id_to_commod_id, paths, pathidx_to_edgelist):
    num_edges = len(paths.keys())
    num_paths = len(set(path_id_to_commod_id.values()))
    
    # set up edge to edge idx, these are original edges in the aggregated graph
    edge_dict = dict()
    for edge_idx, edge in enumerate(paths.keys()):
        edge_dict[edge] = edge_idx

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

def r2_lp():
    pass



def reconciliation_lp(r2_solution_dict, u_meta, v_meta, G, agg_to_orig_nodes, orig_to_agg_node):
    reconciliation_outfile = 'reconciliation_out.txt'
    os.remove(reconciliation_outfile)

    # create a subgraph to represent flows from u_meta to v_meta
    G_u_meta_v_meta = nx.DiGraph()
    nodes_in_u_meta, nodes_in_v_meta = set(), set()

    v_hat_ins, v_hat_outs = v_hat_dict(G_agg)
    u_hat_in, v_hat_out = v_hat_ins[u_meta], v_hat_outs[v_meta]

    all_u_meta_v_meta_inter_edges = [(u,v,G[u][v][CAPACITY]) for u in agg_to_orig_nodes[u_meta] \
                                     for v in G.successors(u) if orig_to_agg_node[v]==v_meta]
    interedge_to_idx = {(u,v): idx for idx, (u, v, _) in enumerate(all_u_meta_v_meta_inter_edges)}

    for (u, v, capacity) in all_u_meta_v_meta_inter_edges:
        if u not in nodes_in_u_meta:
            nodes_in_u_meta.add(u)
            G_u_meta_v_meta.add_node(u)
        if v not in nodes_in_v_meta(v):
            nodes_in_v_meta.add(v)
            G_u_meta_v_meta.add_node(v)
        G_u_meta_v_meta.add_edge(u, v, capacity)

    r2_u_meta_sol_dict, r2_v_meta_sol_dict = r2_solution_dict[u_meta], r2_solution_dict[v_meta]

    # only use meta commodity solution from R2 for flows from u_meta to v_meta
    shared_meta_commodities = set(r2_u_meta_sol_dict.keys()).intersection(set(r2_v_meta_sol_dict.keys()))
    shared_meta_commodities = [commod for commod in shared_meta_commodities if sum([v==v_hat_out for (u,v), _ in r2_u_meta_sol_dict[commod]]) > 0]
    shared_meta_commodities = [commod for commod in shared_meta_commodities if sum([u==u_hat_in for (u,v), _ in r2_v_meta_sol_dict[commod]]) > 0]

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

        # sum outflow from u <= f_{k_local}
        u_outflows_dict = defaultdict()
        for (u,v), outflow in u_meta_flows:
            if v==v_hat_out:
                u_outflows_dict[u] += outflow
        u_outflows_dict = {u: outflow if outflow > 0.0 else 0.0 for u, outflow in u_outflows_dict.items()}

        for u in nodes_in_u_meta:
            outflow = u_outflows_dict[u]
            total_u_meta_flows[k_meta] += outflow

            outgoing_inter_edges_to_idx = [interedge_to_idx[(u,v)] for v in G_u_meta_v_meta.successors(u)]
            print(f"Adding outflow constraint: (sum of commods for edges {outgoing_inter_edges_to_idx}) <= {outflow}")
            m.addConstr(quicksum(commodity_vars[edge, k_local] for edge in outgoing_inter_edges_to_idx) <= outflow, f"outflow [u {u}]][k_meta {k_meta}]")

        # sum of inflow to v <= f_{k_local}
        v_inflows_dict = defaultdict()
        for (u,v), inflow in v_meta_flows:
            if u==u_hat_in: 
                v_inflows_dict[v] += inflow
        v_inflows_dict = {v: inflow if inflow > 0.0 else 0.0 for v, inflow in v_inflows_dict.items()}

        for v in nodes_in_v_meta:
            inflow = v_inflows_dict[v]
            total_v_meta_flows[k_meta] += inflow

            incoming_inter_edges_to_idx = [interedge_to_idx[(u,v)] for u in G_u_meta_v_meta.predecessors(v)]
            print(f"Adding inflow constraint: (sum of commods for edges {incoming_inter_edges_to_idx}) <= {inflow}")
            m.addConstr(quicksum(commodity_vars[edge, k_local] for edge in incoming_inter_edges_to_idx) <= inflow, f"inflow [v {v}]][k_meta {k_meta}]")

    # objective is to maximize flow
    obj = quicksum(commodity_vars)
    m.setObjective(obj, GRB.MAXIMIZE)

    return LpSolver(m, None, reconciliation_outfile), G_u_meta_v_meta, shared_meta_commodities, nodes_in_u_meta, nodes_in_v_meta

def r3_lp():
    pass

if __name__ == '__main__':
    #G = toy_network_1()
    G = toy_network_2()
    tm = generate_uniform_tm(G)
    num_clusters = int(np.sqrt(len(G.nodes)))
    iter_id = 0

    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict, hash_for_clusterid = construct_subproblems(G, tm, num_clusters=num_clusters)

    #print("G clusters dict", [(k, G_clusters_dict[k].nodes()) for k in G_clusters_dict])
    print("agg_edge_dict", agg_edge_dict)
    
    # bundle capacities on inter_edges
    edge_to_bundlecap = bundle_cap(G, agg_edge_dict)
    print("edge_to_bundlecap", edge_to_bundlecap)

    # select paths for r1, this iteration
    paths = path_meta(G, G_agg, num_clusters, edge_to_bundlecap, 0)
    print("paths", paths)
    
    r1_solver, r1_path_to_commod, pathidx_to_edgelist, commodidx_to_info = r1_lp(G, paths, agg_commodities_dict, edge_to_bundlecap)
    print(r1_solver.solve_lp(Method.BARRIER))
    print(r1_solver._model.objVal)
    #print(get_solution_as_mat(r1_solver._model, r1_path_to_commod, paths, pathidx_to_edgelist))
    print("solution as dict", get_solution_as_dict(r1_solver._model, pathidx_to_edgelist, commodidx_to_info, r1_path_to_commod))
