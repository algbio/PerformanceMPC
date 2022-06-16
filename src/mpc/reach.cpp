#include <mpc/graph.h>
#include <mpc/naive.h>
#include <mpc/reach.h>
#include <memory>
#include <vector>

bool naive_reach(int a, int b, Graph &g) {
	std::vector<bool> visited(g.n+1);
	auto dfs = [&g, &b, &visited](auto dfs, int s)->bool {
		if(visited[s])
			return false;
		visited[s] = 1;
		if(s == b)
			return true;
		for(auto &u:g.edge_out[s]) {
			if(dfs(dfs, u))
				return true;
		}
		return false;
	};
	return dfs(dfs, a);
}
std::unique_ptr<reachability_idx> graph_reachability(Graph &g) {
	auto mf = pflowk2(g);
	auto pc = minflow_reduction_path_recover_faster(*mf);
	return std::make_unique<reachability_idx>(g, pc);
}
