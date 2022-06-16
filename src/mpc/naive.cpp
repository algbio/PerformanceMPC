#include "naive.h"
#include "graph.h"
#include <cassert>
#include <algorithm>
#include <array>
#include <vector>
#include <queue>
#include <stack>
#include <functional>
#include <limits>
#include <utility>
#include <iostream>
#include <memory>

// Input graph should have valid and satisfied minflow
void minflow_maxflow_reduction(Flowgraph<Edge::Minflow> &fg, std::function<void(Flowgraph<Edge::Maxflow>&)> maxflow_solver) {
	assert(is_valid_minflow(fg));
	Flowgraph<Edge::Maxflow> fg_red(fg.n, fg.source, fg.sink);
	struct reduction_edge {
		Edge::Maxflow &e_reduction;
		Edge::Minflow &e_original;
		bool reverse;
	};
	std::vector<reduction_edge> v;
	int flow = 0;
	for(auto &[u,e]:fg.edge_out[fg.source])
		flow += e->flow;
	for(int i=1; i<=fg.n; i++) {
		for(auto &[u,e]:fg.edge_out[i]) {
			if(e->flow > e->demand) {
				auto *e2 = fg_red.add_edge(i, u);
				e2->capacity = e->flow - e->demand;
				v.push_back({*e2, *e, false});
			}
		}
	}
	for(int i=1; i<=fg.n; i++) {
		for(auto &[u,e]:fg.edge_out[i]) {
			auto e2 = fg_red.add_edge(u, i);
			e2->capacity = flow;
			v.push_back({*e2, *e, true});
		}
	}
	maxflow_solver(fg_red);
	for(auto &u:v) {
		if(u.reverse) {
			u.e_original.flow += u.e_reduction.flow;
		} else {
			u.e_original.flow -= u.e_reduction.flow;
		}
	}
}

// Dinitz’ Algorithm: The Original Version and Even’s Version 233
// Implementation of DA by Cherkassky
void maxflow_solve_edmonds_karp_DMOD(Flowgraph<Edge::Maxflow> &fg) {
	while(true) {
		std::vector<int> vis(fg.n+1, 0);
		std::vector<int> dist(fg.n+1, std::numeric_limits<int>::max());
		std::queue<int> q;
		q.push(fg.sink);
		dist[fg.sink] = 0;
		vis[fg.sink] = 1;
		while(!q.empty()) {
			int cur = q.front();
			q.pop();
			for(auto &[u,e]:fg.edge_out[cur]) {
				if(e->flow == 0 || vis[u])
					continue;
				vis[u] = 1;
				dist[u] = dist[cur]+1;
				q.push(u);
			}
			for(auto &[u,e]:fg.edge_in[cur]) {
				if(e->flow == e->capacity || vis[u])
					continue;
				vis[u] = 1;
				dist[u] = dist[cur]+1;
				q.push(u);
			}
		}
		if(dist[fg.source] == std::numeric_limits<int>::max()) {
			break;
		}
		std::fill(vis.begin(), vis.end(), 0);
		struct pe {
			Edge::Maxflow *e;
			bool reverse;
		};
		std::vector<pe> path;
		auto dfs = [&vis, &fg, &dist, &path](auto dfs, int s)->bool {
			if(vis[s])
				return false;
			vis[s] = 1;
			for(auto &[u,e]:fg.edge_out[s]) {
				if(e->capacity <= e->flow || dist[s]-1 != dist[u])
					continue;
				path.push_back({e, false});
				if(dfs(dfs, u))
					return true;
				path.pop_back();
			}
			for(auto &[u,e]:fg.edge_in[s]) {
				if(e->flow == 0 || dist[s]-1 != dist[u])
					continue;
				path.push_back({e, true});
				if(dfs(dfs, u))
					return true;
				path.pop_back();
			}
			if(s == fg.sink) {
				int e = std::numeric_limits<int>::max();
				for(auto u:path) {
					if(u.reverse)
						e = std::min(u.e->flow, e);
					else
						e = std::min(u.e->capacity-u.e->flow, e);
				}
				for(auto u:path) {
					if(u.reverse)
						u.e->flow -= e;
					else
						u.e->flow += e;
				}
				return true;
			}
			return false;
		};
		while(dfs(dfs, fg.source)) {
			std::fill(vis.begin(), vis.end(), 0);
			path.clear();
		}
	}
}

void maxflow_solve_edmonds_karp(Flowgraph<Edge::Maxflow> &fg) {
	while(true) {
		std::vector<std::pair<int,Edge::Maxflow*>> visited(fg.n+1);
		std::queue<int> q;
		q.push(fg.source);
		while(!q.empty()) {
			int cur = q.front();
			q.pop();
			if(cur == fg.sink) {
				break;
			}
			for(auto &[u,e]:fg.edge_out[cur]) {
				if(e->flow >= e->capacity || visited[u].first)
					continue;
				visited[u].first = cur;
				visited[u].second = e;
				q.push(u);
			}
			for(auto &[u,e]:fg.edge_in[cur]) {
				if(e->flow == 0 || visited[u].first)
					continue;
				visited[u].first = -cur;
				visited[u].second = e;
				q.push(u);
			}
		}
		if(!visited[fg.sink].first)
			break;
		int cur = fg.sink;
		int delta_flow = std::numeric_limits<int>::max();
		while(cur != fg.source) {
			if(visited[cur].first > 0) {
				delta_flow = std::min(delta_flow, visited[cur].second->capacity-visited[cur].second->flow);
			} else {
				delta_flow = std::min(delta_flow, visited[cur].second->flow);
			}
			cur = abs(visited[cur].first);
		}
		cur = fg.sink;
		while(cur != fg.source) {
			if(visited[cur].first > 0) {
				visited[cur].second->flow += delta_flow;
			} else {
				visited[cur].second->flow -= delta_flow;
			}
			cur = abs(visited[cur].first);
		}
	}
}

// O(pathlen+|E|)
path_cover minflow_reduction_path_recover_faster(Flowgraph<Edge::Minflow> &fg) {
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original graph
	std::vector<int> visited(fg.n+1);
	std::vector<std::vector<int>> cover;
	std::vector<std::vector<std::pair<int, Edge::Minflow*>>::iterator> edge_ptr(fg.n+1);
	for(int i=1; i<=fg.n; i++)
		edge_ptr[i] = fg.edge_out[i].begin();
	while(true) {
		std::fill(visited.begin(), visited.end(), 0);
		std::vector<int> path;
		auto dfs = [&fg, &visited, &path, &edge_ptr](auto dfs, int s)->bool {
			if(s == fg.sink) 
				return true;
			visited[s] = 1;
			while(edge_ptr[s] < fg.edge_out[s].end()) {
				auto &[u,e] = *(edge_ptr[s]);
				if(visited[u] || e->flow == 0) {
					edge_ptr[s]++;
					continue;
				}
				if(dfs(dfs, u)) {
					e->flow--;
					path.push_back(u);
					return true;
				}
			}
			return false;
		};
		if(!dfs(dfs, fg.source))
			break;
		path.push_back(fg.source);
		std::vector<int> real_path;
		for(int i=1; i<path.size()-1; i+=2) {
			real_path.push_back(v_r(path[i]));
		}
		std::reverse(real_path.begin(), real_path.end());
		cover.push_back(real_path);
	}

	for(int i=1; i<=fg.n; i++) {
		for(auto &[u, e]:fg.edge_out[i])
			assert(e->flow == 0);
	}
	return cover;
}

path_cover minflow_reduction_path_recover_fast(Flowgraph<Edge::Minflow> &fg) {
	std::vector<int> visited(fg.n+1);
	std::vector<std::stack<std::vector<int>*>> stk(fg.n+1);
	auto dfs = [&stk, &fg, &visited](auto dfs, int s) {
		if(visited[s])
			return;
		visited[s] = 1;
		for(auto &[u, e]:fg.edge_out[s]) {
			dfs(dfs, u);
		}
		for(auto &[u, e]:fg.edge_in[s]) {
			while(e->flow) {
				e->flow--;
				if(s == fg.sink) {
					auto v = new std::vector<int>();
					v->push_back(s);
					stk[s].push(v);
				}
				stk[s].top()->push_back(u);
				stk[u].push(stk[s].top());
				stk[s].pop();
			}
		}
	};
	dfs(dfs, fg.source);
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original graph
	path_cover cover;
	while(!stk[fg.source].empty()) {
		auto path = stk[fg.source].top();
		stk[fg.source].pop();
		std::vector<int> real_path;
		for(int i=1; i<path->size()-1; i+=2) {
			real_path.push_back(v_r((*path)[i]));
		}
		std::reverse(real_path.begin(), real_path.end());
		cover.push_back(real_path);
		delete path;
	}
	for(int i=1; i<=fg.n; i++) {
		for(auto &[u, e]:fg.edge_out[i])
			assert(e->flow == 0);
	}
	return cover;
}

path_cover minflow_reduction_path_recover(Flowgraph<Edge::Minflow> &fg) {
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original graph
	std::vector<int> visited(fg.n+1);
	std::vector<std::vector<int>> cover;
	while(true) {
		std::fill(visited.begin(), visited.end(), 0);
		std::vector<int> path;
		auto dfs = [&fg, &visited, &path](auto dfs, int s)->bool {
			if(s == fg.sink) 
				return true;
			visited[s] = 1;
			for(auto &[u, e]:fg.edge_out[s]) {
				if(visited[u] || e->flow == 0)
					continue;
				if(dfs(dfs, u)) {
					e->flow--;
					path.push_back(u);
					return true;
				}
			}
			return false;
		};
		if(!dfs(dfs, fg.source))
			break;
		path.push_back(fg.source);
		std::vector<int> real_path;
		for(int i=1; i<path.size()-1; i+=2) {
			real_path.push_back(v_r(path[i]));
		}
		std::reverse(real_path.begin(), real_path.end());
		cover.push_back(real_path);
	}

	for(int i=1; i<=fg.n; i++) {
		for(auto &[u, e]:fg.edge_out[i])
			assert(e->flow == 0);
	}
	return cover;
}

// Find augmenting paths 1 by 1 from residual graph with dfs
void naive_minflow_solve(Flowgraph<Edge::Minflow> &fg) {
	std::vector<bool> visited(fg.n+1);
	while(true) {
		std::fill(visited.begin(), visited.end(), 0);
		auto dfs = [&fg, &visited](auto dfs, int s)->bool {
			if(s == fg.sink)
				return true;
			visited[s] = 1;
			for(auto &[u, e]:fg.edge_out[s]) {
				if(visited[u] || e->demand >= e->flow)
					continue;
				if(dfs(dfs, u)) {
					e->flow--;
					return true;
				}
			}
			for(auto &[u, e]:fg.edge_in[s]) {
				if(visited[u])
					continue;
				if(dfs(dfs, u)) {
					e->flow++;
					return true;
				}
			}
			return false;
		};
		if(!dfs(dfs, fg.source))
			break;
	}
}

std::unique_ptr<Flowgraph<Edge::Minflow>> greedy_minflow_reduction_sparsified(Graph &g, std::function<int(int)> node_weight) {
	struct ggraph {
		int n;
		struct edge {
			int flow = 0;
		};
		std::vector<std::vector<std::pair<int, edge>>> edge_in;
		ggraph(int n):n(n) {
			edge_in.resize(n+1);
		}
	};
	ggraph gg(g.n);
	std::vector<int> topo;
	{std::vector<bool> visited(g.n+1);
	auto dfs = [&gg, &g, &visited, &topo](auto dfs, int s) {
		if(visited[s])
			return;
		visited[s] = 1;
		for(auto &u:g.edge_out[s]) {
			dfs(dfs, u);
			gg.edge_in[u].push_back({s,{0}});
		}
		topo.push_back(s);
	};
	for(int i=1; i<=g.n; i++)
		dfs(dfs, i);
	}

	struct Node_flow {
		int source;
		int sink;
		int flow;
	};
	std::vector<Node_flow> node_flow(g.n+1);
	std::vector<int> max_len(g.n+1);
	std::vector<bool> not_covered(g.n+1, 1);
	while(true) {
		std::vector<std::pair<int, ggraph::edge*>> to(g.n+1);
		std::fill(max_len.begin(), max_len.end(), 0);
		std::pair<int, int> best_node = {0,0};
		std::vector<bool> v2(g.n+1);
		for(auto s:topo) {
			v2[s] = 1;
			if(not_covered[s])
				max_len[s]++;
			if(max_len[s] > best_node.second)
				best_node = {s, max_len[s]};
			if(max_len[s] == 0) { 
				   continue;
			}
			for(auto &[u,e]:gg.edge_in[s]) { 
				if(max_len[s] > max_len[u]) {
					max_len[u] = max_len[s];
					to[u] = {s, &e};
				}
			}
		};
		if(best_node.first == 0)
			break;
		int cur = best_node.first;
		node_flow[cur].source++;
		std::vector<bool> reach(g.n+1);
		for(auto &[u, e]:gg.edge_in[cur]) {
			reach[u] = 1;
		}
		while(to[cur].first != 0) {
			// Sparsify
			not_covered[cur] = 0;
			to[cur].second->flow++;
			node_flow[cur].flow++;
			cur = to[cur].first;
			std::vector<std::pair<int, ggraph::edge>> sparsified;
			for(auto &[u, e]:gg.edge_in[cur]) {
				if(e.flow > 0 || !reach[u])
					sparsified.push_back({u,{e.flow}});
				reach[u] = 1;
			}
			gg.edge_in[cur].clear();
			for(auto u:sparsified)
				gg.edge_in[cur].push_back(u);
		}
		not_covered[cur] = 0;
		node_flow[cur].flow++;
		node_flow[cur].sink++;
	}
	// Reduce to minflow
	int source = g.n*2+1;
	int sink = g.n*2+2;
	auto fgo = std::make_unique<Flowgraph<Edge::Minflow>>(g.n*2+2, source, sink);
	auto v_in = [](int v){return v*2-1;};
	auto v_out = [](int v){return v*2;};
	for(int i=1; i<=g.n; i++) {
		for(auto &[u,e]:gg.edge_in[i]) {
			auto *e2 = fgo->add_edge(v_out(u), v_in(i));
			e2->flow = e.flow;
		}
		Edge::Minflow *e = fgo->add_edge(v_in(i), v_out(i));
		e->demand = 1;
		e->flow = node_flow[i].flow;
		e = fgo->add_edge(source, v_in(i));
		e->demand = 0;
		e->flow = node_flow[i].source;
		e = fgo->add_edge(v_out(i), sink);
		e->demand = 0;
		e->flow = node_flow[i].sink;
	}
	assert(is_valid_minflow(*fgo));
	return fgo;
}

std::unique_ptr<Flowgraph<Edge::Minflow>> greedy_minflow_reduction(Graph &g, std::function<int(int)> node_weight) {
	Flowgraph<Edge::Minflow> tfg = {g.n, 0, 0};
	std::vector<int> topo;
	std::vector<bool> visited(g.n+1);
	auto dfs = [&tfg, &g, &visited, &topo](auto dfs, int s) {
		if(visited[s])
			return;
		visited[s] = 1;
		std::vector<int> qq;
		for(auto &u:g.edge_out[s]) {
			qq.push_back(u);
		}
		for(auto &u:qq) {
			dfs(dfs, u);
			tfg.add_edge(s, u);
		}
		topo.push_back(s);
	};
	std::vector<int> qq;
	for(int i=1; i<=g.n; i++) {
		qq.push_back(i);
	}
	for(auto u:qq) {
		dfs(dfs, u);
	}
	std::reverse(topo.begin(), topo.end());

	struct Node_flow {
		int source;
		int sink;
		int flow;
	};
	std::vector<Node_flow> node_flow(g.n+1);
	std::vector<int> max_len(g.n+1);
	std::vector<std::pair<int, Edge::Minflow*>> from(g.n+1);
	auto &not_covered = visited;
	while(true) {
		std::fill(max_len.begin(), max_len.end(), 0);
		std::pair<int, int> best_node = {0,0};
		for(auto s:topo) {
			if(not_covered[s])
				max_len[s]++;
			if(max_len[s] > best_node.second)
				best_node = {s, max_len[s]};
			if(max_len[s] == 0) {
				continue;
			}
			for(auto &[u,e]:tfg.edge_out[s]) {
				if(max_len[s] > max_len[u]) {
					max_len[u] = max_len[s];
					from[u] = {s, e};
				}
			}
		}
		if(best_node.first == 0)
			break;
		int cur = best_node.first;
		node_flow[cur].sink++;
		while(from[cur].first != 0) {
			not_covered[cur] = 0;
			from[cur].second->flow++;
			node_flow[cur].flow++;
			cur = from[cur].first;
		}
		not_covered[cur] = 0;
		node_flow[cur].flow++;
		node_flow[cur].source++;
	}
	// Reduce to minflow
	int source = g.n*2+1;
	int sink = g.n*2+2;
	auto fgo = std::make_unique<Flowgraph<Edge::Minflow>>(g.n*2+2, source, sink);
	auto v_in = [](int v){return v*2-1;};
	auto v_out = [](int v){return v*2;};
	for(int i=1; i<=g.n; i++) {
		for(auto &[u,e]:tfg.edge_out[i]) {
			auto *e2 = fgo->add_edge(v_out(i), v_in(u));
			e2->flow = e->flow;
		}
		Edge::Minflow *e = fgo->add_edge(v_in(i), v_out(i));
		e->demand = 1;
		e->flow = node_flow[i].flow;
		e = fgo->add_edge(source, v_in(i));
		e->demand = 0;
		e->flow = node_flow[i].source;
		e = fgo->add_edge(v_out(i), sink);
		e->demand = 0;
		e->flow = node_flow[i].sink;
	}
	assert(is_valid_minflow(*fgo));
	return fgo;
}

std::unique_ptr<Flowgraph<Edge::Minflow>> naive_minflow_reduction(Graph &g, std::function<int(int)> node_weight) {
	int source = g.n*2+1;
	int sink = g.n*2+2;
	auto fgo = std::make_unique<Flowgraph<Edge::Minflow>>(g.n*2+2, source, sink);
	auto v_in = [](int v){return v*2-1;};
	auto v_out = [](int v){return v*2;};
	for(int i=1; i<=g.n; i++) {
		for(auto &u:g.edge_out[i]) {
			fgo->add_edge(v_out(i), v_in(u));
		}
		Edge::Minflow *e = fgo->add_edge(v_in(i), v_out(i));
		e->demand = node_weight(i);
		e->flow = node_weight(i);
		e = fgo->add_edge(source, v_in(i));
		e->demand = 0;
		e->flow = node_weight(i);
		e = fgo->add_edge(v_out(i), sink);
		e->demand = 0;
		e->flow = node_weight(i);
	}
	return fgo;
}

bool is_valid_minflow(Flowgraph<Edge::Minflow> &fg) {
	for(int i=1; i<=fg.n; i++) {
		int total_out = 0;
		for(auto &[u,e]:fg.edge_out[i]) {
			if(e->flow < e->demand) {
				std::cout << "Demand not satisfied" << std::endl;
				return false;
			}
			total_out += e->flow;
		}
		int total_in = 0;
		for(auto &[u,e]:fg.edge_in[i])
			total_in += e->flow;
		if(i != fg.sink && i != fg.source && total_in != total_out) {
				std::cout << "Flow conservation not satisfied " << i << " " << total_in << "/" << total_out<< std::endl;
			return false;
		}
	}
	return true;
}

bool is_valid_cover(std::vector<std::vector<int>> &cover, Graph &g) {
	std::vector<int> visited(g.n+1);
	for(auto path:cover) {
		for(auto &u:path) {
			if(u < 1 || u > g.n) {
				std::cout << " oob " << std::endl;
				return false;
			}
			visited[u] = 1;
		}
		for(int i=1; i<path.size(); i++)
			if(!g.has_edge(path[i-1], path[i])) {
				std::cout << " no edge " << std::endl;
				return false;
			}
	}
	for(int i=1; i<=g.n; i++)
		if(!visited[i]) {
				std::cout << " no visit all "<< i << std::endl;
			return false;
		}
	return true;
}
