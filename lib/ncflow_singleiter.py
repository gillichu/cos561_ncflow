from LpSolver import *
import gurobipy
from gurobipy import GRB, Model, quicksum
from itertools import tee

def path_to_edge_list(path):
	a, b = tee(path)
	next(b, None)
	return zip(a, b)


def solve_lp():
	pass 



def r1_lp(G_meta, paths_dict, agg_commodity_dict):
	r1_outfile = 'r1_out.txt'
	commodities = []
	edge_to_paths = defaultdict(list)
	r1_path_to_commodities = dict()
	r1_commodity_to_pathids = defaultdict(list)
	
	path_idx = 0

	for key in agg_commodity_dict.keys():
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
		num_paths.append(len(path_ids))
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
	for u, v, c_e in G_meta.edges.data('cap'):
		if (u, v) in edge_to_paths: 
			paths = edge_to_paths[(u, v)]
			constr_vars = [path_vars[p] for p in paths]
			m.addConstr(quicksum(constr_vars) <= c_e)

	return LpSolver(m, None, r1_outfile) 

def r2_lp():
	pass

def reconciliation_lp():
	pass

def r3_lp():
	pass
