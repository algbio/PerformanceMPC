#include "graph.h"
#include <utility>
#include <memory>
#include <iostream>

std::unique_ptr<Graph> random_dag(int n, int m, int seed) {
	auto g = std::make_unique<Graph>(n);
	std::vector<int> topo(n);
	for(int i=1; i<=n; i++)
		topo[i-1] = i;
    std::mt19937 rg(seed);
	std::shuffle(topo.begin(), topo.end(), rg);

	int edges = 0;
	std::uniform_int_distribution<int> dst(1, n);
	while(edges < m) {
		int a = dst(rg);
		int b = dst(rg);
		if(topo[a-1] >= topo[b-1] || g->has_edge(a, b))
			continue;
		g->add_edge(a, b);
		edges++;
	}
	return g;
}

std::unique_ptr<Graph> random_x_partite(int x, int n, int m, int seed) {
	auto g = std::make_unique<Graph>(n*x);
    std::mt19937 rg(seed);
	int edges = 0;
	std::uniform_int_distribution<int> dst1(1, n);
	for(int i=0; i<x-1; i++) {
		int ofs = i*n;
		std::uniform_int_distribution<int> dst2(ofs+n+1, x*n);
		for(int j=0; j<std::min(m, (x-1-i)*n); j++) {
			int a = ofs+dst1(rg);
			int b = dst2(rg);
			if(g->has_edge(a, b))
				continue;
			g->add_edge(a, b);
		}
	}
	return g;
}

// Complete binary tree directed from root
std::unique_ptr<Graph> binary_tree(int depth, bool reverse) {
	int n = (1<<depth);
	auto g = std::make_unique<Graph>(n);
	for(int i=1; i<=n; i++) {
		if(reverse) {
			if(i*2 <= n)
				g->add_edge(i*2, i);
			if(i*2+1 <= n)
				g->add_edge(i*2+1, i);
		} else {
			if(i*2 <= n)
				g->add_edge(i, i*2);
			if(i*2+1 <= n)
				g->add_edge(i, i*2+1);
		}
	}
	return g;
}

// N nodes connected to "tunnel" of N nodes and then again N nodes (width N, pathlen = N*N)
std::unique_ptr<Graph> funnel_gen(int n) {
	auto g = std::make_unique<Graph>(n*3); 
	for(int i=1; i<=n; i++) {
		g->add_edge(i, n+1);
		if(i < n)
			g->add_edge(n+i, n+i+1);
		g->add_edge(n+n, n+n+i);
	}
	return g;
}

std::unique_ptr<Graph> random_x_chain(int x, int n, int m, int seed) {
    std::mt19937 rg(seed);
											   
	auto g = random_dag(n, m, seed);
	std::vector<int> v;
	for(int i=1; i<=n; i++)
		v.push_back(i);
	std::shuffle(v.begin(), v.end(), rg);
	std::vector<int> topo;
	std::vector<bool> visited(g->n+1);
	auto dfs = [&g, &rg, &topo, &visited](int s, auto dfs) -> void {
		if(visited[s])
			return;
		visited[s] = 1;
		std::vector<int> v(g->edge_out[s]);
		std::shuffle(v.begin(), v.end(), rg);
		for(auto u:v) {
			dfs(u, dfs);
		}
		topo.push_back(s);
	};
	for(auto u:v) {
		dfs(u, dfs);
	}
	std::reverse(topo.begin(), topo.end());
	std::vector<int> last(x+1);
	std::uniform_int_distribution<> dst(1, x);
	for(int i=0; i<n; i++) {
		int lol = dst(rg);
		if(last[lol] && !g->has_edge(last[lol], topo[i])) {
			g->add_edge(last[lol], topo[i]);
		}
		last[lol] = topo[i];
	}
	return g;
}
