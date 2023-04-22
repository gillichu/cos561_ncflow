import re
from LpSolver import *
from gurobipy import GRB, Model, quicksum
from path import *
from create_subproblems import *
from itertools import tee, combinations

def nodelist_to_edge_list(path):
    a, b = tee(path)
    next(b, None)
    return zip(a, b)


def solve_lp():
    pass 


def r1_lp(G, paths_dict, agg_commodities_dict):
    ### expects input:
    ### paths_dict[(u_meta, v_meta)] = [([e1, e2, e3, ..., eN], mincap)]
    print("path dictionary", paths_dict)

    r1_outfile = 'r1_out.txt'
    commodities = []
    r1_path_to_commodities = dict()
    r1_commodity_to_pathids = defaultdict(list)
    edge_to_pathids = defaultdict(list)
    cap_list = nx.get_edge_attributes(G, 'capacity')
    print("cap list", cap_list)
    
    path_idx = 0

    for key in agg_commodities_dict.keys():
        commodity_idx, pathinfo = key
        s_k, t_k, d_k = pathinfo
        print("looking up:", s_k, t_k)

        # get single path between each pair of meta nodes
        # should only have 1, since selected inter-cluster edges
        path_nodelist = paths_dict[(s_k, t_k)][0] 
        # original node edges
        for edge in nodelist_to_edge_list(path_nodelist):
            print("original edge", edge, "demand:", d_k, "capacity", cap_list[edge])
            edge_to_pathids[edge].append(path_idx)
            
        # get the path to commodity (for r3)
        r1_path_to_commodities[path_idx] = commodity_idx
        r1_commodity_to_pathids[commodity_idx].append(path_idx)
        
        commodities.append((commodity_idx, d_k, [path_idx]))
        path_idx += 1
    
    m = Model('max-flow: R1')
    
    # create a variable for each path
    path_variables = m.addVars(path_idx, vtype=GRB.CONTINUOUS, lb=0.0, name='f')

    # set objective
    obj = quicksum(path_variables)
    m.setObjective(obj, GRB.MAXIMIZE)

    # add demand constraints
    for _, d_k, path_ids in commodities:
        # sum of all path variables for commodity k (only one) should be <= commodity k's demand (d_k)
        m.addConstr(quicksum(path_variables[p] for p in path_ids) <= d_k)

    # add edge capacity constraints 
    for edge in edge_to_pathids.keys():

        # get all paths on this edge
        path_indices = edge_to_pathids[edge]

        # get edge capacity
        c_e = cap_list[edge]

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
        for edge in nodelist_to_edge_list(all_paths[p]):
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

    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict, hash_for_clusterid = construct_subproblems(G, tm, num_clusters=num_clusters)

    print("G clusters dict", [(k, G_clusters_dict[k].nodes()) for k in G_clusters_dict])

    # select paths for r1, this iteration
    paths = path_meta(G, G_agg, num_clusters, agg_edge_dict, 0)
    
    r1_solver, r1_path_to_commod = r1_lp(G, paths, agg_commodities_dict)
    print(r1_solver.solve_lp(Method.BARRIER))
    print(r1_solver._model.objVal)
    #print(get_solver_results(r1_solver._model, G_agg, r1_path_to_commod, all_r1_paths))
