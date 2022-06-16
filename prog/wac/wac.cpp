#include <mpc/lemon.h>
#include <mpc/naive.h>
#include <mpc/graph.h>
#include <mpc/pflow.h>
#include <mpc/utils.cpp>
#include <mpc/antichain.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <nlohmann/json.hpp>
#include <signal.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <algorithm>
#include <limits>
#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

std::vector<std::pair<std::function<void(Flowgraph<Edge::Minflow>&)>, std::string>> mifsol = {
	{lemon_ns, "lemon_ns"},
	{lemon_ns_minlen, "lemon_ns_minlen"}, 
	{lemon_cs, "lemon_cs"},
	{lemon_cs_minlen, "lemon_cs_minlen"},
	{lemon_cc, "lemon_cc"},
	{lemon_cc_minlen, "lemon_cc_minlen"},
	{lemon_caps, "lemon_caps"},
	{lemon_caps_minlen, "lemon_caps_minlen"},
	{naive_minflow_solve, "naive_minflow_solve"}};
std::vector<std::pair<std::function<void(Flowgraph<Edge::Maxflow>&)>, std::string>> mafsol =
	{{lemon_preflow, "lemon_preflow"},
	{maxflow_solve_edmonds_karp, "maxflow_solve_edmonds_karp"},
	{maxflow_solve_edmonds_karp_DMOD, "maxflow_solve_edmonds_karp_DMOD"}};


void graph_info(nlohmann::json &j, Graph &g) {
	j["graph"]["n"] = g.n;
	int edges = 0;
	for(int i=1; i<=g.n; i++)
		edges += g.edge_out[i].size();
	j["graph"]["m"] = edges;
}

std::pair<Graph*, std::vector<int>> read_graph2(std::string filename) {
	std::ifstream input(filename);
	if(!input.good()) {
		std::cerr << "Failed to open file: " << filename << "\n";
		exit(1);
	}
	int n, m;
	input >> n >> m;
	Graph *g = new Graph(n);
	std::vector<int> weight(n+1);
	for(int i=1; i<=n; i++) {
		input >> weight[i];
	}
	for(int i=0; i<m; i++) {
		int a, b;
		input >> a >> b;
		g->add_edge(a, b);
	}
	return {g, weight};
}

int main(int argc, char** args) {
	struct rlimit rlm;
	int ret = getrlimit(RLIMIT_STACK, &rlm);
	if(ret != 0)
		exit(1);
	if(rlm.rlim_cur != RLIM_INFINITY) {
		std::cerr << "Warning: RLIMIT_STACK not set to infinity, but: " << rlm.rlim_cur << ", e.g. ulimit -s unlimited" << std::endl;
	}

	CLI::App app{R"(Input graph file format:
	<n m> number of nodes and edges
	n lines <w_i> weight of node i
	m lines <a b> edge from a to b

Output format:
	<s> chain size
	s lines, each containing vertex of the chain

Solver names correspond to
	naive_minflow_solve = Find decrementing paths with dfs
Implementations from LEMON graph library:
	ns = NetworkSimplex
	cs = CostScaling
	cc = CycleCanceling
	caps = CapacityScaling
	minlen = Set cost of every edge to 1 in Mincostflow reduction
From LEMON documentation:
	"In general, NetworkSimplex and CostScaling are the most efficient implementations. NetworkSimplex is usually the fastest on relatively small graphs (up to several thousands of nodes) and on dense graphs, while CostScaling is typically more efficient on large graphs (e.g. hundreds of thousands of nodes or above), especially if they are sparse. However, other algorithms could be faster in special cases. For example, if the total supply and/or capacities are rather small, CapacityScaling is usually the fastest algorithm (without effective scaling)."
			)"};

	std::string filepath;
	std::string solver;
	std::string reduction;
	app.add_option("-f,--file", filepath, "Path to input file")->mandatory();
	std::string solvers_list = "";
	for(auto &u:mifsol) {
		if(!solvers_list.empty())
			solvers_list.append(", ");
		solvers_list.append(u.second);
	}
	app.add_option("-s,--solver", solver, "Solver to use {"+solvers_list+"}")->mandatory();
	app.add_option("-r,--reduction", reduction, "Reduction to use {naive}")->mandatory();
	CLI11_PARSE(app, argc, args);

	if(reduction != "naive") {
		std::cerr << "Unknown reduction: " << reduction << "\n";
		exit(1);
	}

	auto solver1 = std::find_if(mifsol.begin(), mifsol.end(), [&solver](auto u){return solver == u.second;});
	//auto solver2 = std::find_if(mafsol.begin(), mafsol.end(), [&solver_s](auto u){return solver_s == u.second;});
	//if(solver1 == mifsol.end() && solver2 == mafsol.end()) {
	if(solver1 == mifsol.end()) {
		std::cerr << "Unknown slover: " << solver << std::endl;
		exit(1);
	}
	auto graph = read_graph2(filepath);

	auto fg = naive_minflow_reduction(*graph.first, [&graph](int i){return graph.second[i];});
	solver1->first(*fg);
	auto chain = maxantichain_from_minflow(*fg);
	std::cout << chain.size() << "\n";
	for(auto u:chain)
		std::cout << u << "\n";
	std::cout << "\n";
}
