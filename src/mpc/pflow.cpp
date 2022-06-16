#include "pflow.h"
#include <iterator>
#include "naive.h"
#include "graph.h"
#include <queue>
#include <vector>
#include <deque>
#include <algorithm>
#include <limits>
#include <iostream>
#include <list>

std::unique_ptr<Flowgraph<Edge::Minflow>> pflowk3(Graph &g) {
	std::vector<int> topo;
	// Topo order O(n+m) OK
	{
		std::vector<bool> visited(g.n+1);
		auto s1 = [&topo, &g, &visited](int s, auto dfs) {
			if(visited[s])
				return;
			visited[s] = 1;
			for(auto u:g.edge_in[s]) {
				dfs(u, dfs);
			}
			topo.push_back(s);
		};
		for(int i=1; i<=g.n; i++)
			s1(i, s1);
	}
	std::vector<int> topo_idx(g.n+1);
	for(int i=0; i<topo.size(); i++) {
		topo_idx[topo[i]] = i;
	}
	auto v_in = [](int v){return v*2-1;};
	auto v_out = [](int v){return v*2;};
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original graph
	auto fgo = std::make_unique<Flowgraph<Edge::Minflow>>(g.n*2+2, g.n*2+1, g.n*2+2);
	auto &fg = *fgo;
	std::vector<int> pathid(g.n+1); // some path id for node x
	std::vector<std::vector<int>> pathids(g.n+1); // path ids of a node x
	int paths = 0; // how many paths in total
	std::vector<int> path_ends; // Nodes in which flowpaths end
	std::vector<int> layer(fg.n+1); // layer of node x
	std::vector<std::vector<int>> layer_v(2); // nodes of layer x
	layer[fg.source] = std::numeric_limits<int>::min();
	layer[fg.sink] = std::numeric_limits<int>::max();
	std::vector<std::pair<int, std::vector<std::pair<int, Edge::Minflow*>>::iterator>> path_s_v(fg.n+1, {-1, {}}); 
	std::vector<bool> visited2(fg.n+2); 
	std::vector<std::vector<int>> new_paths(g.n+1); 
	std::vector<std::pair<int, Edge::Minflow*>> visited(fg.n+1); 
	for(int i=0; i<topo.size(); i++) {
		int vi = topo[i];
		auto *edge = fg.add_edge(v_in(vi), v_out(vi));
		edge->demand = 1;
		edge->flow = 1;
		auto *edge2 = fg.add_edge(fg.source, v_in(vi));
		edge2->flow = 1;
		edge = fg.add_edge(v_out(vi), fg.sink);
		edge->flow = 1;
		// Sparsify
		std::vector<int> sparsify(paths, -1);
		for(auto u:g.edge_in[vi]) {
			sparsify[pathid[u]] = std::max(sparsify[pathid[u]], topo_idx[u]);
		}
		for(auto u:sparsify) {
			if(u == -1)
				continue;
			fg.add_edge(v_out(topo[u]), v_in(vi));
		}
		// Find a dec path
		std::vector<std::queue<int>> q(layer_v.size());
		q[q.size()-1].push(v_in(vi));
		std::vector<int> visited_v;
		visited[fg.source] = {1, nullptr};
		visited[v_in(vi)] = {fg.source, edge2};
		std::vector<int> visited_reset;
		visited_reset.push_back(v_in(vi));
		layer[v_in(vi)] = std::numeric_limits<int>::max();
		layer[v_out(vi)] = std::numeric_limits<int>::max();
		bool path_found = false;
		int lowest_lr = std::numeric_limits<int>::max();
		for(int g=q.size()-1; g>=0; g--) {
			while(!q[g].empty() && !path_found) {
				int s = q[g].front();
				q[g].pop();
				if(s == fg.sink)
					continue;
				lowest_lr = std::min(lowest_lr, layer[s]);
				visited_v.push_back(s);
				for(auto u:fg.edge_out[s]) {
					if(visited[u.first].first || u.second->flow <= u.second->demand)
						continue;
					visited[u.first] = {s, u.second};
					visited_reset.push_back(u.first);
					if(u.first == fg.sink) {
						// Dec path found
						path_found = true;
						break;
					}
					q[layer[u.first]].push(u.first);
				}
				for(auto u:fg.edge_in[s]) {
					if(visited[u.first].first)
						continue;
					visited[u.first] = {-s, u.second};
					visited_reset.push_back(u.first);
					q[layer[u.first]].push(u.first);
				}
			}
		}
		if(!path_found) {
			pathids[vi].push_back(paths);
			pathid[vi] = paths;
			paths++;
			path_ends.push_back(v_out(vi));
		}
		if(path_found) {
			// fix flow
			int cur = fg.sink;
			while(cur != fg.source) {
				auto *e = visited[cur].second;
				cur = visited[cur].first;
				if(cur > 0) {
					e->flow--;
				} else {
					e->flow++;
				}
				cur = abs(cur);
			}
		}
		for(auto u:visited_reset)
			visited[u] = {0, nullptr};
		// We can assume that a path was found
		// Update path ids
		std::vector<int> used_path(paths);
		if(path_found){
			if(!path_found)
				std::cout << lowest_lr << std::endl;
			std::vector<int> path_s_v_reset;
			auto s2 = [&fg, &layer, &v_out, &vi,&v_r, &path_ends, &path_s_v, &path_s_v_reset, &new_paths, &pathid, &pathids](int s, int pid, auto fs) {
				if(s == fg.sink)
					return;
				if(s%2 == 0) {
					// only add paths if is v_out node
					pathid[v_r(s)] = pid;
					new_paths[v_r(s)].push_back(pid);
				}
				if(path_s_v[s].first == -1) {
					path_s_v[s].second = fg.edge_out[s].begin();
					path_s_v[s].first = path_s_v[s].second->second->flow;
					path_s_v_reset.push_back(s);
				}
				while(path_s_v[s].first == 0) {
					path_s_v[s].second++;
					path_s_v[s].first = path_s_v[s].second->second->flow;
				}
				path_s_v[s].first--;
				if(path_s_v[s].second->first == fg.sink) {
					path_ends[pid] = s;
				}
				fs(path_s_v[s].second->first, pid, fs);
			};
			std::vector<int> visited2_reset;
			auto s1 = [&fg, &s2, &new_paths, &used_path, &v_r, &v_out, &path_ends, &pathids, &visited2, &visited2_reset, &layer, &lowest_lr](int s, auto fs)  {
				if(visited2[s] || layer[s] < lowest_lr || s == fg.sink || s == fg.source)
					return;
				visited2[s] = 1;
				visited2_reset.push_back(s);
				for(auto &u:fg.edge_in[s]) {
					fs(u.first, fs);
				}
				// s2 is called in topological order
				for(auto &u:pathids[v_r(s)]) {
					if(used_path[u])
						continue;
					used_path[u] = 1;
					s2(s, u, s2);
				}
				if(s%2==0) {
					pathids[v_r(s)] = new_paths[v_r(s)];
					new_paths[v_r(s)].clear();
				}
			};
			s1(v_in(vi), s1);
			s1(v_out(vi), s1);
			for(int g=lowest_lr; g<layer_v.size(); g++) {
				for(auto &u:layer_v[g]) {
					s1(u, s1);
				}
			}
			for(auto u:path_s_v_reset)
				path_s_v[u] = {-1, {}};
			for(auto u:visited2_reset)
				visited2[u] = 0;
		}
		// Update layers
		if(lowest_lr == std::numeric_limits<int>::max()) {
			// No other vertex visited
			lowest_lr = 0;
			if((unsigned long) lowest_lr+1 == layer_v.size()) {
				// top lr visited create new layer for v_out
				layer_v.push_back({});
			}
			layer[v_in(vi)] = lowest_lr;
			layer_v[lowest_lr].push_back(v_in(vi));
			layer[v_out(vi)] = lowest_lr+1;
			layer_v[lowest_lr+1].push_back(v_out(vi));
			continue;
		}
		if((unsigned long) lowest_lr+1 == layer_v.size()) {
			// top lr visited create new layer for v_out
			layer_v.push_back({});
		}
		for(auto u:visited_v) {
			if(layer[u] != lowest_lr)
				layer[u] = -1;
		}
		for(int g=lowest_lr+1; g<layer_v.size(); g++) {
			std::vector<int> v;
			for(auto u:layer_v[g]) {
				if(layer[u] == g)
					v.push_back(u);
			}
			layer_v[g] = v;
		}
		for(auto u:visited_v) {
			if(v_r(u) == vi)
				continue;
			if(layer[u] != lowest_lr) {
				layer_v[lowest_lr].push_back(u);
				layer[u] = lowest_lr;
			}
		}
		layer[v_in(vi)] = lowest_lr;
		layer_v[lowest_lr].push_back(v_in(vi));
		layer[v_out(vi)] = lowest_lr+1;
		layer_v[lowest_lr+1].push_back(v_out(vi));
		// Merge if no path ends at l
		bool should_merge = true;
		for(auto u:path_ends) {
			if(layer[u] == lowest_lr) {
				should_merge = false;
			}
		}
		if(should_merge && lowest_lr > 0) {
			// merge l+1 to l ...
			for(int g=lowest_lr; g<layer_v.size(); g++) {
				for(auto u:layer_v[g]) {
					layer[u] = g-1;
					layer_v[g-1].push_back(u);
				}
				layer_v[g].clear();
			}
			layer_v.pop_back();
		}
	}
	return fgo;
}

std::unique_ptr<Flowgraph<Edge::Minflow>> pflowk2(Graph &g) {
	std::vector<int> topo;
	topo.reserve(g.n);
	// Topo order O(n+m) OK
	{
		std::vector<bool> visited(g.n+1);
		auto s1 = [&topo, &g, &visited](int s, auto dfs) {
			if(visited[s])
				return;
			visited[s] = 1;
			for(auto &u:g.edge_in[s]) {
				dfs(u, dfs);
			}
			topo.push_back(s);
		};
		for(int i=1; i<=g.n; i++)
			s1(i, s1);
	}
	std::vector<int> topo_idx(g.n+1);
	for(int i=0; i<topo.size(); i++) {
		topo_idx[topo[i]] = i;
	}
	auto v_in = [](int v){return v*2-1;};
	auto v_out = [](int v){return v*2;};
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original graph
	auto fgo = std::make_unique<Flowgraph<Edge::Minflow>>(g.n*2+2, g.n*2+1, g.n*2+2);
	auto &fg = *fgo;
	std::vector<int> backlink(g.n+1);
	std::vector<int> backlink_id(g.n+1);
	int backlink_cnt = 0;
	std::vector<int> newlink(g.n+1);
	std::vector<int> backlink_first(g.n+1);
	std::vector<int> layer(fg.n+1); // layer of node x
	std::vector<std::vector<int>> layer_v(2);
	layer[fg.source] = std::numeric_limits<int>::min();
	layer[fg.sink] = std::numeric_limits<int>::max();
	std::vector<std::pair<int, std::vector<std::pair<int, Edge::Minflow*>>::iterator>> path_s_v(fg.n+1, {-1, {}});
	std::vector<bool> visited2(fg.n+2); 
	std::vector<std::vector<int>> new_paths(g.n+1);
	std::vector<std::pair<int, Edge::Minflow*>> visited(fg.n+1); 
	for(int i=0; i<topo.size(); i++) {
		const int vi = topo[i];
		auto *edge = fg.add_edge(v_in(vi), v_out(vi));
		edge->demand = 1;
		edge->flow = 1;
		edge = fg.add_edge(v_out(vi), fg.sink);
		edge->flow = 1;
		edge = fg.add_edge(fg.source, v_in(vi));
		edge->flow = 1;
		// Sparsify
		std::vector<int> sparsify(backlink_cnt, -1);
		for(auto &u:g.edge_in[vi]) {
			auto blid = backlink_id[backlink[u]];
			sparsify[blid] = std::max(sparsify[blid], topo_idx[u]);
		}
		for(auto u:sparsify) {
			if(u == -1)
				continue;
			fg.add_edge(v_out(topo[u]), v_in(vi));
		}
		// Find a dec path
		std::vector<std::queue<int>> q(layer_v.size());
		q[q.size()-1].push(v_in(vi));
		std::vector<int> visited_v;
		visited[fg.source] = {1, nullptr};
		visited[v_in(vi)] = {fg.source, edge};
		std::vector<int> visited_reset;
		visited_reset.push_back(v_in(vi));
		layer[v_in(vi)] = std::numeric_limits<int>::max();
		layer[v_out(vi)] = std::numeric_limits<int>::max();
		bool path_found = false;
		int lowest_lr = std::numeric_limits<int>::max();
		for(int g=q.size()-1; g>=0; g--) {
			while(!q[g].empty() && !path_found) {
				int s = q[g].front();
				q[g].pop();
				lowest_lr = g;
				visited_v.push_back(s);
				for(auto &u:fg.edge_out[s]) {
					if(visited[u.first].first || u.second->flow <= u.second->demand)
						continue;
					visited[u.first] = {s, u.second};
					visited_reset.push_back(u.first);
					if(u.first == fg.sink) {
						// Dec path found
						path_found = true;
						break;
					}
					q[layer[u.first]].push(u.first);
				}
				for(auto &u:fg.edge_in[s]) {
					if(visited[u.first].first)
						continue;
					visited[u.first] = {-s, u.second};
					visited_reset.push_back(u.first);
					q[layer[u.first]].push(u.first);
				}
			}
		}
		if(path_found) {
			int cur = fg.sink;
			while(cur != fg.source) {
				auto *e = visited[cur].second;
				cur = visited[cur].first;
				if(cur > 0) {
					e->flow--;
				} else {
					e->flow++;
				}
				cur = abs(cur);
			}
		}
		for(auto u:visited_reset)
			visited[u] = {0, nullptr};
		if(!path_found) {
			backlink[vi] = vi;
			backlink_id[vi] = backlink_cnt++;
			backlink_first[vi] = 1;
		}
		if(lowest_lr == std::numeric_limits<int>::max()) {
			lowest_lr = 0;
			if((unsigned long) lowest_lr+1 == layer_v.size()) {
				layer_v.push_back({});
			}
			layer[v_in(vi)] = lowest_lr;
			layer_v[lowest_lr].push_back(v_in(vi));
			layer[v_out(vi)] = lowest_lr+1;
			layer_v[lowest_lr+1].push_back(v_out(vi));
			continue;
		}
		if((unsigned long) lowest_lr+1 == layer_v.size()) {
			layer_v.push_back({});
		}
		for(auto u:visited_v) {
			if(layer[u] != lowest_lr)
				layer[u] = -1;
		}
		for(int g=lowest_lr+1; g<layer_v.size(); g++) {
			std::vector<int> v;
			for(auto u:layer_v[g]) {
				if(layer[u] == g)
					v.push_back(u);
			}
			layer_v[g] = v;
		}
		for(auto u:visited_v) {
			if(v_r(u) == vi)
				continue;
			if(layer[u] != lowest_lr)
				layer_v[lowest_lr].push_back(u);
			layer[u] = lowest_lr;
		}
		layer[v_in(vi)] = lowest_lr;
		layer_v[lowest_lr].push_back(v_in(vi));
		layer[v_out(vi)] = lowest_lr+1;
		layer_v[lowest_lr+1].push_back(v_out(vi));
		bool should_merge = true;
		{
			std::vector<int> path_s_v_reset;
			std::vector<std::pair<int,int>> vp1;
			std::vector<std::pair<int,int>> vp2;
			auto s2 = [&backlink, &newlink,&vi, &vp1, &vp2, &lowest_lr,&backlink_first, &fg, &layer, &v_out,&v_in, &v_r, &path_s_v, &path_s_v_reset](int s, int bli, auto fs)->int {
				if(s == fg.sink)
					return -1;
				if(layer[s] > lowest_lr) {
					if(v_r(s) != vi){
						vp2.push_back({backlink[v_r(s)], v_r(s)});
					}
					return v_r(s);
				}
				if(path_s_v[s].first == -1) {
					path_s_v[s].second = fg.edge_out[s].begin();
					assert(path_s_v[s].second != fg.edge_out[s].end());
					path_s_v[s].first = path_s_v[s].second->second->flow;
					path_s_v_reset.push_back(s);
				}
				while(path_s_v[s].first == 0) {
					path_s_v[s].second++;
					assert(path_s_v[s].second != fg.edge_out[s].end());
					path_s_v[s].first = path_s_v[s].second->second->flow;
				}
				path_s_v[s].first--;
				int nwl = fs(path_s_v[s].second->first, bli, fs);
				if(v_r(s) != bli && v_in(v_r(s)) == s)
					vp1.push_back({v_r(s), bli});
				if(v_r(s) != bli && v_in(v_r(s)) == s) 
					newlink[v_r(s)] = nwl;
				return nwl;
			};
			std::vector<int> visited2_reset;
			auto s1 = [&fg, &s2,&should_merge, &lowest_lr,&backlink_first, &v_in,  &v_r, &v_out, &visited2, &visited2_reset, &layer](int s, auto fs)  {
				if(visited2[s] || layer[s] < lowest_lr || s == fg.sink || s == fg.source)
					return;
				visited2[s] = 1;
				visited2_reset.push_back(s);
				if(layer[s] == lowest_lr && lowest_lr > 0) {
					for(auto &u:fg.edge_out[s]) {
						if(u.first == fg.sink && u.second->flow > 0) {
							should_merge = false;
						}
					}
				}
				for(auto &u:fg.edge_in[s]) {
					fs(u.first, fs);
				}
				if((s == v_out(v_r(s)) && layer[v_in(v_r(s))] < lowest_lr) || (s == v_in(v_r(s)) && backlink_first[v_r(s)]) ) {
					s2(s, v_r(s), s2);
				}
			};
			for(auto u:layer_v[lowest_lr]) {
				s1(u, s1);
			}
			for(auto u:path_s_v_reset)
				path_s_v[u] = {-1, {}};
			for(auto u:visited2_reset)
				visited2[u] = 0;
			for(auto u:vp1) {
				backlink[u.first] = u.second; 
			}
			for(auto u:vp2) {
				newlink[u.first] = u.second;
			}
			for(auto u:path_s_v_reset)
				path_s_v[u] = {-1, {}};
			for(auto u:visited2_reset)
				visited2[u] = 0;
			for(int g=lowest_lr; g<layer_v.size(); g++) {
				for(auto u:layer_v[g]) {
					if(g > lowest_lr) {
						if(layer[v_out(backlink[v_r(u)])] != layer[v_in(v_r(u))] && !backlink_first[v_r(u)]) {
							backlink[v_r(u)] = newlink[backlink[v_r(u)]];
						}
					}
					if(u == v_out(v_r(u)) && layer[v_in(v_r(u))] != layer[v_out(v_r(u))]) {
						backlink_id[v_r(u)] = backlink_id[backlink[v_r(u)]];
					}
				}
			}
		}
		if(should_merge && lowest_lr > 0) {
			for(int g=lowest_lr; g<layer_v.size(); g++) {
				std::vector<std::pair<int,int>> qq;
				if(g == lowest_lr) {
					for(auto u:layer_v[g]) {
						if(u == v_in(v_r(u)) && layer[v_in(backlink[v_r(u)])] == lowest_lr-1) {
							qq.push_back({v_r(u), backlink[backlink[v_r(u)]]});
						}
					}
					for(auto u:qq)
						backlink[u.first] = u.second;
				}
				for(auto u:layer_v[g]) {
					layer[u] = g-1;
					layer_v[g-1].push_back(u);
				}
				layer_v[g].clear();
			}
			layer_v.pop_back();
		}
	}
	return fgo;
}
