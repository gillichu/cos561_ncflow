from traffic_matrix import *
from wrapper import solver_wrapper
from benchmarks import *
from preprocess import * 

def time_series(G, tm,num_tms, rel_delta_abs_mean, rel_delta_std, num_cluster,max_num_iters,outfile,test_pf4=False):
    tms = generate_sequence_of_tms(tm,num_tms, rel_delta_abs_mean, rel_delta_std)
    demand_results, time_results = [],[]
    if test_pf4: pf4_demands, pf4_times = [],[]
    num_nodes = len(G.nodes())
    for i in range(0,len(tms)):
        tm = tms[i]
        print(i)
        print(tm)
        G, tm, num_nodes = preprocess_tz_files("", G, tm, num_nodes)
        total_demand, time_taken = solver_wrapper(G, tm, max_num_iters, int(num_cluster))
        demand_results.append(total_demand)
        time_results.append(time_taken)
        if test_pf4: 
            pf4_demand, pf4_time = pf4_run(G, tm)
            pf4_demands.append(pf4_demand)
            pf4_times.append(pf4_time)

    outfile = 'time_series_' + outfile + '_result.txt'
    with open(outfile, 'w+') as f:
        f.write(str(num_cluster) + "\n")
        f.write(str(demand_results) + "\n")
        if test_pf4: f.write("PF4: "+ str(pf4_demands) + "\n")
        f.write(str(time_results) + "\n")
        if test_pf4: f.write("PF4: "+ str(pf4_times) + "\n")
    return demand_results,time_results

if __name__ == '__main__':
    max_num_iter = 10
    num_tms = 3
    perturb_mean = 1
    perturb_std = 0.5
    G = read_graphml('/Users/cloverzsq/Desktop/cos561_ncflow/topologies/Uninett2010.graphml')
    tm = generate_uniform_tm(G)
    num_cluster = int(np.sqrt(len(G.nodes)))
    demands, times = time_series(G,tm,num_tms,perturb_mean,perturb_std,num_cluster,max_num_iter,"test1",True)
    print(demands)
    print(times)
    # num_nodes = len(G.nodes())
    # G, tm, num_nodes = preprocess_tz_files("", G, tm, num_nodes)
    # total_demand, time_taken = solver_wrapper(G, tm, max_num_iter, int(num_cluster))
    # print(total_demand)
    # print(time_taken)
