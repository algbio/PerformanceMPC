#include <vector>
#include <stack>
#include <mpc/naive.h>
#include <mpc/graph.h>
#include <mpc/antichain.h>


antichain maxantichain_from_minflow(Flowgraph<Edge::Minflow> &mf) {
	std::vector<int> visited(mf.n+1);
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original graph
	antichain mac;
	auto dfs = [&mac, &v_r, &mf, &visited](auto &dfs, int s) {
		if(visited[s])
			return;
		assert(s != mf.sink);
		visited[s] = 1;
		for(auto &[u,e] : mf.edge_out[s]) {
			if(e->flow > e->demand) {
				dfs(dfs, u);
			}
		}
		for(auto &[u, e] : mf.edge_in[s]) {
			dfs(dfs, u);
		}
	};
	dfs(dfs, mf.source);
	auto dfs2 = [&mac, &v_r, &mf, &visited](auto &dfs, int s) {
		if(visited[s] != 1)
			return;
		visited[s] = 2;
		for(auto &[u,e] : mf.edge_out[s]) {
			if(e->flow > e->demand) {
				dfs(dfs, u);
			}
			if(e->flow == e->demand && e->demand >= 1 && !visited[u]) {
				mac.push_back(v_r(s));
				visited[u] = 3; 
			}

		}
		for(auto &[u, e] : mf.edge_in[s]) {
			dfs(dfs, u);
		}
	};
	dfs2(dfs2, mf.source);
	return mac;
}

bool is_antichain(antichain &ac, Graph &g) {
	std::vector<bool> visited(g.n+1), antichain(g.n+1);
	for(auto u:ac)
		antichain[u] = 1;
	auto dfs = [&g, &visited, &antichain](auto &dfs, int s)->bool {
		if(visited[s])
			return false;
		if(antichain[s]) {
			return true;
		}
		visited[s] = 1;
		for(auto u:g.edge_out[s]) {
			if(dfs(dfs, u))
				return true;
		}
		return false;
	};
	for(auto u:ac) {
		antichain[u] = 0;
		if(dfs(dfs, u))
			return false;
		antichain[u] = 1;
		for(int i=1; i<=g.n; i++)
			visited[i] = 0;
	}
	return true;
}
