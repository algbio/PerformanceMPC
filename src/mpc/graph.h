#pragma once
#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <iostream>
#include <memory>
#include <iostream>
#include <deque>

struct Graph {

	std::vector<std::vector<int>> edge_in, edge_out;
	int n;

	Graph(int n):n(n) {
		edge_in.resize(n+1);
		edge_out.resize(n+1);
	}
	Graph(const Graph &g):n(g.n) {
		edge_in.resize(n+1);
		edge_out.resize(n+1);
		for(int i=1; i<=n; i++)  {
			edge_in[i] = g.edge_in[i];
			edge_out[i] = g.edge_out[i];
		}
	}

	int count_edges() {
		int edges = 0;
		for(int i=1; i<=n; i++)
			edges += edge_out[i].size();
		return edges;
	}

	void add_edge(int a, int b) {
		edge_out[a].push_back(b);
		edge_in[b].push_back(a);
	}

	bool has_edge(int a, int b) {
		for(auto &u : edge_out[a])
			if(u == b)
				return true;
		return false;
	}


	void print() {
		int m=0;
		for(int i=1; i<=n; i++) {
			m += edge_out[i].size();
		}
		std::cout << n << " " << m << "\n";
		for(int i=1; i<=n; i++) {
			for(auto u:edge_out[i]) {
				std::cout << i << " " << u << "\n";
			}
			m += edge_out[i].size();
		}
		std::cout << std::endl;
	}

	void delete_edge(int a, int b) {
		edge_out[a].erase(std::find(edge_out[a].begin(), edge_out[a].end(), b));
		edge_in[b].erase(std::find(edge_in[b].begin(), edge_in[b].end(), a));
	}

};

namespace Edge {
	struct Minflow {
		int flow = 0;
		int demand = 0;
	};
	struct Maxflow {
		int flow = 0;
		int capacity = 0;
	};
	struct Maxflow_skew : Maxflow {
		int flow = 0;
		int capacity = 0;
		Maxflow_skew *reverse;
	};
};

template<typename EdgeT> struct Flowgraph {


	int n, sink, source;
	std::vector<std::vector<std::pair<int, EdgeT*>>> edge_in, edge_out;
	std::deque<EdgeT> edges;

	Flowgraph(int n, int source, int sink):n(n),sink(sink),source(source) {
		edge_in.resize(n+1);
		edge_out.resize(n+1);
	}


	EdgeT* add_edge(int a, int b) {
		edges.push_back({});
		EdgeT *e = &*edges.rbegin();
		edge_out[a].push_back({b, e});
		edge_in[b].push_back({a, e});
		return e;
	}

	Flowgraph(const Flowgraph &fg):n(fg.n),sink(fg.sink),source(fg.source) {
		edge_in.resize(n+1);
		edge_out.resize(n+1);
		for(int i=1; i<=n; i++)
			for(auto &u:fg.edge_out[i]) {
				auto *e = add_edge(i, u.first);
				e->demand = u.second->demand;
				e->flow = u.second->flow;
			}
	}

	bool has_edge(int a, int b) {
		for(auto &u : edge_out[a])
			if(u.first == b)
				return true;
		return false;
	}

	EdgeT* find_edge(int a, int b) {
		auto is_b = [&b](std::pair<int, EdgeT*> x) {return x.first==b;};
		return std::find_if(edge_out[a].begin(), edge_out[a].end(), is_b)->second;
	}

	void delete_edge(int a, int b) {
		auto is_b = [&b](std::pair<int, EdgeT*> x) {return x.first==b;};
		auto e1 = std::find_if(edge_out[a].begin(), edge_out[a].end(), is_b);
		delete e1->second;
		edge_out[a].erase(e1);
		auto is_a = [&a](std::pair<int, EdgeT*> x) {return x.first==a;};
		edge_in[b].erase(std::find_if(edge_in[b].begin(), edge_in[b].end(), is_a));
	}

	~Flowgraph() {
	}
};

std::unique_ptr<Graph> funnel_gen(int n);
std::unique_ptr<Graph> random_dag(int n, int m, int seed);
std::unique_ptr<Graph> random_x_partite(int x, int n, int m, int seed);
std::unique_ptr<Graph> random_x_chain(int k, int n, int m, int seed);
std::unique_ptr<Graph> binary_tree(int depth, bool reverse);
