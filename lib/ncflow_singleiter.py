import os
import re
from LpSolver import *
from gurobipy import GRB, Model, quicksum
from path import *
from create_subproblems import *
from itertools import tee, combinations
# from r2 import *

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

def r2_lp(current_meta_node, paths_dict, G, G_meta, all_v_hat_in, all_v_hat_out, intra_commods, partition_vector, r1_sol_dict, meta_commodity_dict, r1_sol_mat, meta_to_virt_dict, virt_to_meta_dict):
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

def extract_r2_solution_dict(model, agg_commodities_dict, multi_commodity_list, path_id_to_commod_id, all_paths, virt_to_meta_dict): 
    commod_id_to_meta_commod_id = {commod_key: meta_commod_key[0] for meta_commod_key, c_l in agg_commodities_dict.items() for commod_key, (_, _, _) in c_l}

    meta_sol_dict_def = defaultdict(list)
    mc_id_to_path_id_to_flow = defaultdict(dict)
    active_meta_commodity_dict = {} 
 

    for srcs, targets, _, commod_list in multi_commodity_list:
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
            meta_sol_dict_def[meta_commod_key] += [
                    (edge, var.x) for edge in path_to_edge_list(all_paths[p])
                ]
 
    return {commod: meta_sol_dict_def[commod] if commod in meta_sol_dict_def else [] for commod in active_meta_commodity_dict.values()}



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

    return LpSolver(m, None, reconciliation_outfile), G_u_meta_v_meta, shared_meta_commodities, nodes_in_u_meta, nodes_in_v_meta, all_u_meta_v_meta_inter_edges

def extract_reconciliation_solution_dict(G, r2_solution_dict, agg_to_orig_nodes, orig_to_agg_node, all_v_hat_in, all_v_hat_out):
    num_clusters = len(agg_to_orig_nodes)
    all_metanode_pairs = [(u_meta, v_meta) for u_meta in range(num_clusters) for v_meta in range(num_clusters) if u_meta != v_meta]

    reconciliation_solutions_dicts = {}
    for (u_meta, v_meta) in all_metanode_pairs:
        reconciliation_solver, G_u_meta_v_meta, shared_meta_commodities, nodes_in_u_meta, nodes_in_v_meta, all_u_meta_v_meta_inter_edges = reconciliation_lp(r2_solution_dict, u_meta, v_meta, G, agg_to_orig_nodes, orig_to_agg_node,meta_to_virt_dict)
        reconciliation_solver.solve_lp(Method.BARRIER)
        
        model = reconciliation_solver._model
        l = []
        for var in model.getVars():
            match = re.match(r'{}\[(\d+),(\d+)\]'.format('f'), var.varName)
            edge_idx, k_local = int(match.group(1)), int(match.group(2))
            edge = list(G_u_meta_v_meta.edges)[edge_idx]#all_u_meta_v_meta_inter_edges[edge_idx]

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

def r3_lp():
    pass

if __name__ == '__main__':
    # G = toy_network_1()
    G = toy_network_2()
    tm = generate_uniform_tm(G)
    num_clusters = int(np.sqrt(len(G.nodes)))
    iter_id = 0

    G_agg, agg_edge_dict, agg_to_orig_nodes, orig_to_agg_node, G_clusters_dict, agg_commodities_dict,clusters_commodities_dict, hash_for_clusterid, meta_to_virt_dict, virt_to_meta_dict = construct_subproblems(G, tm, num_clusters=num_clusters)

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
    r2_solution_dict = {}

    for current_meta_node in G_agg.nodes:
        paths_dict = path_r2(current_meta_node,G_clusters_dict[current_meta_node],all_v_hat_in, all_v_hat_out, agg_to_orig_nodes)
        intra_commods = clusters_commodities_dict[current_meta_node]
        v_hat_in, v_hat_out = all_v_hat_in[current_meta_node], all_v_hat_out[current_meta_node]
        r2_solver, multi_commodity_list, path_id_to_commod_id, all_paths = r2_lp(current_meta_node, paths_dict, G, G_agg, v_hat_in, v_hat_out, intra_commods, orig_to_agg_node, r1_sol_dict, agg_commodities_dict, r1_sol_mat, meta_to_virt_dict, virt_to_meta_dict)
        r2_solver.solve_lp(Method.BARRIER)

        model = r2_solver._model
        print('current_meta_node ',current_meta_node)
        r2_solution_dict[current_meta_node] = extract_r2_solution_dict(model, agg_commodities_dict, multi_commodity_list, path_id_to_commod_id, all_paths, virt_to_meta_dict)


    # format: (u_meta, v_meta) --> (k_meta, (s_k, t_k, d_k)) --> (true edge, flow)
    reconciliation_solutions_dicts = extract_reconciliation_solution_dict(G, r2_solution_dict, agg_to_orig_nodes, orig_to_agg_node,  all_v_hat_in, all_v_hat_out)
    print('Reconciliation', reconciliation_solutions_dicts)

    print(r1_sol_dict)

    # reconciliation_meta_out_flows = defaultdict(lambda: defaultdict(float))
    # reconciliation_meta_in_flows = defaultdict(lambda: defaultdict(float))

    # # meta_commod_id -> meta_node_id -> out/in flow
    # for (u_meta, v_meta), reconciliation_sol_dict in reconciliation_solutions_dicts.items():
    #     # k_meta to flow
    #     after_recon_flow = dict()
    #     for meta_commod, flow_list in reconciliation_sol_dict.items():
    #         total_flow = 0.0
    #         for (u, v), l in flow_list:
    #             if orig_to_agg_node[u] == u_meta and orig_to_agg_node[v] == v_meta:
    #                 total_flow += l
    #         k_meta = meta_commod[0]
    #         reconciliation_meta_out_flows[k_meta][u_meta] += total_flow
    #         reconciliation_meta_in_flows[k_meta][v_meta] += total_flow
    #         after_recon_flow[k_meta] = total_flow
    #     print('(u_meta, v_meta)',(u_meta, v_meta))
    #     print('after_recon_flow',after_recon_flow)

    # print('reconciliation_meta_out_flows',reconciliation_meta_out_flows)
    # print('reconciliation_meta_in_flows',reconciliation_meta_in_flows)

    