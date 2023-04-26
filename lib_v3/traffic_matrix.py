import numpy as np
from collections import defaultdict
import networkx as nx


def generate_uniform_tm(G, max_demand=10, seed=0):
    '''
    Generate uniform traffic matrix for network G
    '''
    np.random.seed(seed)
    num_nodes = len(G.nodes)
    tm = np.random.rand(num_nodes, num_nodes) * max_demand
    tm = tm.astype(np.float32)
    np.fill_diagonal(tm, 0.0)
    
    return tm

def generate_gravity_tm(G, total_demand, seed=0, sample=False):
    '''
    Generate gravity traffic matrix for network G
    '''
    np.random.seed(seed)
    num_nodes = len(G.nodes)
    tm = np.zeros((num_nodes, num_nodes), dtype=np.float32)

    sccs = nx.strongly_connected_components(self.problem.G)
    for scc in sccs:
        in_cap_sum, out_cap_sum = defaultdict(float), defaultdict(float)
        for u in scc:
            for v in G.predecessors(u):
                in_cap_sum[u] += G[v][u]['capacity']
            for v in G.successors(u):
                out_cap_sum[u] += G[u][v]['capacity']
        in_cap_sum, out_cap_sum = dict(in_cap_sum), dict(out_cap_sum)

        in_total_cap = sum(in_cap_sum.values())
        out_total_cap = sum(out_cap_sum.values())

        for u in scc:
            norm_u = out_cap_sum[u] / out_total_cap
            for v in scc:
                if u == v:
                    continue
                frac = norm_u * in_cap_sum[v] / \
                    (in_total_cap - in_cap_sum[u])
                if sample:
                    tm[u, v] = max(np.random.normal(frac, frac / 4), 0.0)
                else:
                    tm[u, v] = frac

    tm *= total_demand
    return tm

def generate_poisson_tm(G, lam, const_factor, decay=1.0, seed=0):
    '''
    Generate poisson traffic matrix for network G
    '''
    np.random.seed(seed)
    num_nodes = len(G.nodes)

    distances = np.zeros((num_nodes, num_nodes), dtype=np.int)
    dist_iter = nx.shortest_path_length(G)
    for src, dist_dict in dist_iter:
        for target, dist in dist_dict.items():
            distances[src, target] = dist

    tm = np.array([[np.random.poisson(lam * (decay**dist)) for dist in row] for row in distances],
                        dtype=np.float32)
    np.fill_diagonal(tm, 0.0)

    tm *= const_factor
    return tm

def generate_gaussian_tm(G, mean, stddev, decay=1.0, seed=0):
    '''
    Generate gaussian traffic matrix for network G
    '''
    np.random.seed(seed)
    num_nodes = len(G.nodes)

    tm = np.random.normal(mean, stddev, (num_nodes, num_nodes))
    tm[tm < 0.0] = 0.0
    np.fill_diagonal(tm, 0.0)
    
    return tm

def generate_bimodal_tm(G, fraction, low_range, high_range, seed=0):
    '''
    Generate bimodal traffic matrix for network G
    '''
    np.random.seed(seed)
    num_nodes = len(G.nodes)

    tm = np.zeros((num_nodes, num_nodes))
    inds = np.random.choice(2, (num_nodes, num_nodes), p=[
                            fraction, 1 - fraction]).astype('bool')
    tm[inds] = np.random.uniform(low_range[0], low_range[1], np.sum(inds))
    tm[~inds] = np.random.uniform(high_range[0], self.high_range[1], np.sum(~inds))
    np.fill_diagonal(tm, 0.0)
        
    return tm

def perturb_matrix(orig_tm, mean, stddev):
        tm = np.copy(orig_tm)
        tm += np.random.normal(mean, stddev, tm.shape)
        np.fill_diagonal(tm, 0.0)
        tm[tm < 0.0] = 0.0

        return tm

# call this function to generate a sequence of traffic matrices for experiments
# (each traffic matrix is generated based on the one from the previous timestep)
def generate_sequence_of_tms(seed_tm, num_tms, rel_delta_abs_mean, rel_delta_std, seed=1): 
    mean_load = np.mean(seed_tm)
   
    delta_mean = rel_delta_abs_mean * mean_load
    delta_std = rel_delta_std * mean_load

    all_tms = [seed_tm]

    np.random.seed(seed)
    while len(all_tms) < num_tms:

        curr_seed_tm = all_tms[-1]
        # Change each demand by delta mean on average, with a spread of delta_std
        perturb_mean = np.random.normal(delta_mean, delta_std)
        perturb_mean *= np.random.choice([-1, 1])
        tm = perturb_matrix(curr_seed_tm, perturb_mean, delta_std)
        all_tms.append(tm)

    return all_tms