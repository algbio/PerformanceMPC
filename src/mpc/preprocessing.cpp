#include "preprocessing.h"
#include "naive.h"
#include <memory>


std::unique_ptr<contract_tree_recovery_data> contract_tree(Graph &g) {
	std::vector<int> topo;
	{	std::vector<bool> visited(g.n+1);
		auto dfs = [&g, &visited, &topo](auto dfs, int s) {
			if(visited[s])
				return;
			visited[s] = 1;
			for(auto &u:g.edge_out[s]) {
				dfs(dfs, u);
			}
			topo.push_back(s);
		};
		for(int i=1; i<=g.n; i++)
			dfs(dfs, i);
		std::reverse(topo.begin(), topo.end());
	}
	std::vector<bool> visited(g.n+1, false);


	auto result = std::make_unique<contract_tree_recovery_data>();
	result->node_recover.push_back({});
	int ncc = 0;


	std::vector<bool> no_include(g.n+1);
	std::vector<std::pair<int,int>> new_edges;
	std::vector<int> edge_added(g.n+1, 0);
	std::vector<int> mapping(g.n+1); // node id in contracted graph
	
	auto new_id = [&mapping, &result, &ncc](int s)->int {
		if(mapping[s] == 0) {
			mapping[s] = ++ncc;
			result->node_recover.push_back({s});
		}
		return mapping[s];
	};
	auto dfs_down = [&new_id, &edge_added, &result, &visited, &g, &new_edges](int s, auto dfs_down, int root, std::vector<int> &pth)->void {
		visited[s] = true;
		bool leaf = true;
		for(auto u:g.edge_out[s]) {
			if(g.edge_in[u].size() == 1) {
				leaf = false;
			}
		}
		for(auto u:g.edge_out[s]) {
			if(g.edge_in[u].size() == 1) {
				pth.push_back(u);
				dfs_down(u, dfs_down, root, pth);
				pth.pop_back();
			} else {
				if(leaf) {
					new_edges.push_back({new_id(s), new_id(u)});
				} else {
					if(edge_added[u] != root) {
						new_edges.push_back({new_id(root), new_id(u)});
						if(pth.size() > 0)
							result->edge_recover[{new_id(root), new_id(u)}] = {pth.begin(), pth.end()};
						edge_added[u] = root;
					}
				}
			}
		}
		if(s == root)
			new_id(root);
		if(s != root && leaf) {
			new_edges.push_back({new_id(root), new_id(s)});
			if(pth.size() > 1) {
				result->node_recover[new_id(s)] = {pth.begin(), pth.end()};
			}
		}
	};

	for(auto u:topo) {
		if(!visited[u]) {
			std::vector<int> pth;
			dfs_down(u, dfs_down, u, pth);
		}
	}


	result->g = {ncc};
	for(auto u:new_edges) {
		result->g.add_edge(u.first, u.second);
	}
	return result;
}

path_cover recover_contract_pathcover(path_cover &pc, contract_tree_recovery_data &rec) {
	path_cover np = {};
	for(auto u:pc) {
		np.push_back({});
		int pv = -1;
		for(auto s:u) {
			if(pv != -1) {
				if(rec.edge_recover.find({pv, s}) != rec.edge_recover.end()) {
					for(auto v:rec.edge_recover[{pv, s}])
						np.rbegin()->push_back(v);
				}
			}
			for(auto v:rec.node_recover[s]) {
				np.rbegin()->push_back(v);
			}
			pv = s;
		}
	}
	return np;
}

path_cover recover_contract_pathcover(path_cover &pc, std::vector<std::vector<int>> rcv) {
	path_cover np = {};
	for(auto u:pc) {
		np.push_back({});
		for(auto s:u) {
			for(auto v:rcv[s]) {
				np.rbegin()->push_back(v);
			}
		}
	}
	return np;
}

// O(n+m) 
std::unique_ptr<Graph> sparsify_dfs_elegant(Graph &gs) {
	std::vector<std::vector<int>> edge_out_topo_order(gs.n+1);
	std::vector<int> topo;
	{
	std::vector<bool> visited(gs.n+1);
	int ans = 0;
	auto dfs = [&gs, &visited, &topo, &edge_out_topo_order](auto dfs, int s) {
		if(visited[s])
			return;
		visited[s] = 1;
		for(auto &u:gs.edge_out[s]) {
			dfs(dfs, u);
		}
		topo.push_back(s);
	};
	for(int i=1; i<=gs.n; i++)
		dfs(dfs, i);
	std::reverse(topo.begin(), topo.end());
	for(auto u:topo) {
		for(auto uu:gs.edge_in[u]) {
			edge_out_topo_order[uu].push_back(u);
		}
	}
	}
	auto g = std::make_unique<Graph>(gs.n);
	std::vector<bool> visited(g->n+1);
	std::vector<int> last_dfs_reach(g->n+1);
	int dfs_reach_cnt = 1;
	auto dfs = [&last_dfs_reach, &dfs_reach_cnt, &visited, &gs, &g, &edge_out_topo_order](int s, auto dfs) {
		if(visited[s])
			return;
		visited[s] = 1;
		int cur_reach = dfs_reach_cnt++;
		for(auto u:edge_out_topo_order[s]) {
			dfs(u, dfs);
			if(last_dfs_reach[u] < cur_reach) {
				g->add_edge(s, u);
				last_dfs_reach[u] = cur_reach;
			}
		}
	};
	for(int i=0; i<gs.n; i++) {
		dfs(topo[i], dfs);
	}
	return g;
}
