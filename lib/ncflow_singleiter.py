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


def r1_lp(G_meta, paths_dict, agg_commodities_dict):
    r1_outfile = 'r1_out.txt'
    commodities = []
    edge_to_paths = defaultdict(list)
    r1_path_to_commodities = dict()
    r1_commodity_to_pathids = defaultdict(list)
    
    path_idx = 0


    for key in agg_commodities_dict.keys():
        commodity_k, pathinfo = key
        s_k, t_k, d_k = pathinfo

        # get single path between each pair of meta nodes
        paths = paths_dict[(s_k, t_k)]
        assert len(paths) == 1
        path = paths[0]
        path_ids = []

        # get all the edges in the path
        for edge in path_to_edge_list(path):
            edge_to_paths[edge].append(path_idx)
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
    for u, v, c_e in G_meta.edges.data('capacity'):
        # print("G_meta", G_meta.edges.data)
        if (u, v) in edge_to_paths: 
            paths = edge_to_paths[(u, v)]
            constr_vars = [path_variables[p] for p in paths]
            m.addConstr(quicksum(constr_vars) <= c_e)

    return LpSolver(m, None, r1_outfile), r1_path_to_commodities

def get_solver_results(model, G, path_id_2_commod_id, all_paths):
    num_edges = len(G.edges)
    num_paths = len(set(path_id_2_commod_id.values()))
    # set up edge to edge idx
    edge_dict = dict()
    for edge_idx, edge in enumerate(G.edges):
        edge_dict[edge] = edge_idx
    print("Edge_dict", edge_dict)

    sol_mat = np.zeros((num_edges, num_paths), dtype=np.float32)
    for var in model.getVars():
        # match var name back to path
        p = int(re.match(r'f\[(\d+)\]', var.varName).group(1))
        commodity_idx = path_id_2_commod_id[p]
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

    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict = construct_subproblems(G, tm, num_clusters=num_clusters)

    paths_dict = path_meta(G, num_clusters, agg_edge_dict, agg_to_orig_nodes)

    # select paths for r1, this iteration
    r1_paths = dict()
    all_r1_paths = []
    
    # create r1_agg_commodities_dict
    for key in agg_commodities_dict.keys():
        commodity_k, pathinfo = key
        s_k, t_k, d_k = pathinfo
        paths = paths_dict[(s_k, t_k)]
        r1_paths[(s_k, t_k)] = [paths[0]]
        all_r1_paths.append(paths[0])
    
    # add edges back into G_meta. I'm currently doing this incorrectly. I think this is also the reason get_solver_results is confused. 
    for a, b in combinations(G_agg.nodes, 2): 
        if a != b:
            add_bi_edge(G_agg, a, b, 2)

    print("G_agg", G_agg)

    r1_solver, r1_path_to_commod = r1_lp(G_agg, r1_paths, agg_commodities_dict)
    print(r1_solver.solve_lp(Method.BARRIER))
    print(r1_solver._model.objVal)
    #print(get_solver_results(r1_solver._model, G_agg, r1_path_to_commod, all_r1_paths))
