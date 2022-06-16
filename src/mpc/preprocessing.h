#pragma once
#include "graph.h"
#include "naive.h"
#include <map>

std::unique_ptr<Graph> sparsify_dfs_elegant(Graph &gs);

struct contract_tree_recovery_data {
	Graph g;
	std::vector<std::vector<int>> node_recover;
	std::map<std::pair<int,int>, std::vector<int>> edge_recover;
	contract_tree_recovery_data() : g(0) {
	}
};
std::unique_ptr<contract_tree_recovery_data> contract_tree(Graph &g);
path_cover recover_contract_pathcover(path_cover &pc, contract_tree_recovery_data &rec);
