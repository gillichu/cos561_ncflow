import re
from LpSolver import *
from gurobipy import GRB, Model, quicksum
from path import *
from create_subproblems import *
from itertools import tee, combinations

def path_to_edge_list(path):
    a, b = tee(path)
    next(b, None)
    return zip(a, b)


def solve_lp():
    pass 


def r1_lp(paths_dict, agg_commodities_dict):
    ### expects input:
    ### paths_dict[(u_meta, v_meta)] = [([e1, e2, e3, ..., eN], mincap)]

    r1_outfile = 'r1_out.txt'
    commodities = []
    r1_path_to_commodities = dict()
    r1_commodity_to_pathids = defaultdict(list)
    edge_to_paths = defaultdict(list)
    
    path_idx = 0

    for key in agg_commodities_dict.keys():
        commodity_k, pathinfo = key
        s_k, t_k, d_k = pathinfo

        # get single path between each pair of meta nodes
        paths = paths_dict[(s_k, t_k)]
        
        path = paths[0]
        for edge in path_to_edge_list(path):
            edge_to_paths[edge].append(path_idx)

        path_ids = []
        path_ids.append(path_idx)

        # get the path to commodity (for r3)
        r1_path_to_commodities[path_idx] = commodity_k
        r1_commodity_to_pathids[commodity_k].append(path_idx)
        
        commodities.append((commodity_k, d_k, path_ids))
        path_idx += 1
    
    m = Model('max-flow: R1')
    
    # create a variable for each path
    path_variables = m.addVars(path_idx, vtype=GRB.CONTINUOUS, lb=0.0, name='f')

    # set objective
    obj = quicksum(path_variables)
    m.setObjective(obj, GRB.MAXIMIZE)

    # add demand constraints
    for _, d_k, path_ids in commodities:
        # all path ids should be <= commodity k's demand (d_k)
        m.addConstr(quicksum(path_variables[p] for p in path_ids) <= d_k)

    # add edge capacity constraints 
    for u, v in paths_dict: # for each intercluster edge in the selected paths dictionary
        # paths_dict[(u_meta, v_meta)] = [([e1, e2, e3, ..., eN], mincap)]
        _, c_e = paths_dict[(u, v)]
        path_indices = edge_to_paths[(u, v)] 

        # ensure that all paths on a given edge meet the edge constraint
        constr_vars = [path_variables[p] for p in path_indices]
        m.addConstr(quicksum(constr_vars) <= c_e)

    return LpSolver(m, None, r1_outfile), r1_path_to_commodities

def get_solver_results(model, G, path_id_to_commod_id, all_paths):
    num_edges = len(G.edges)
    num_paths = len(set(path_id_to_commod_id.values()))
    # set up edge to edge idx, these are original edges in the aggregated graph
    edge_dict = dict()
    for edge_idx, edge in enumerate(G.edges):
        edge_dict[edge] = edge_idx
    print("Edge_dict", edge_dict)

    sol_mat = np.zeros((num_edges, num_paths), dtype=np.float32)
    for var in model.getVars():
        # match var name back to path
        p = int(re.match(r'f\[(\d+)\]', var.varName).group(1))
        commodity_idx = path_id_to_commod_id[p]
        for edge in path_to_edge_list(all_paths[p]):
            sol_mat[edge_dict[edge], commodity_idx] += var.x
    return sol_mat

def r2_lp():
    pass

def reconciliation_lp():
    pass

def r3_lp():
    pass

if __name__ == '__main__':
    G = toy_network_2()
    tm = generate_uniform_tm(G)
    num_clusters = int(np.sqrt(len(G.nodes)))
    iter_id = 0

    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict = construct_subproblems(G, tm, num_clusters=num_clusters)

    # original edges on the collapsed graph, outputs dict (metanode, metanode) = set of paths
    paths_dict = path_meta(G, G_agg, num_clusters, agg_edge_dict, iter_id)

    # select paths for r1, this iteration
    selected_inter_edges = select_inter_edge(G, agg_edge_dict, iter_id)
    print('selected_inter_edges', selected_inter_edges)
    
    r1_solver, r1_path_to_commod = r1_lp(selected_inter_edges, agg_commodities_dict)
    print(r1_solver.solve_lp(Method.BARRIER))
    print(r1_solver._model.objVal)
    #print(get_solver_results(r1_solver._model, G_agg, r1_path_to_commod, all_r1_paths))
