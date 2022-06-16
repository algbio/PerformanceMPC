#include <algorithm>
#include <memory>
#include <mpc/cc.cpp>
#include <mpc/lemon.h>
#include <mpc/naive.h>
#include <mpc/graph.h>
#include <mpc/pflow.h>
#include <mpc/utils.h>
#include <mpc/transitive.h>
#include <mpc/preprocessing.h>
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
#include "CLI/Validators.hpp"

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

void run_one(Graph &g, unsigned int timeout_sec, unsigned long mem_limit_bytes, std::string reduction_s, std::string solver_s, bool sparsify_dfs_f, bool contract_trees, std::string output_path, std::string cover_decomposition, nlohmann::json &j) {
	j["reduction"]["name"] = reduction_s;
	j["solver"]["name"] = solver_s;
	j["mem_limit_bytes"] = mem_limit_bytes;
	j["time_limit_sec"] = timeout_sec;
	graph_info(j, g);

	stopwatch sw_c = {sw_child};
	auto pid = fork();
	if(pid < 0) {
		exit(-1);
	} else if(pid == 0) {
		rlimit rlm = {mem_limit_bytes, mem_limit_bytes};
		setrlimit(RLIMIT_AS, &rlm);
		signal(SIGALRM, [](int signum) {exit(124);});
		alarm(timeout_sec);
		stopwatch sw_s = {sw_self};
		///////
		j["preprocess"]["sparsify_dfs"]["enabled"] = sparsify_dfs_f;
		j["preprocess"]["contract_tree"]["enabled"] = contract_trees;
		if(sparsify_dfs_f) {
			int cn = g.count_edges();
			j["preprocess"]["sparsify_dfs"]["before"]["m"] = g.count_edges();
			stopwatch sparsify_time = {};
			g = *sparsify_dfs_elegant(g);
			log_time(sparsify_time.total(), j["sparsify_dfs"]["time"]);
			j["preprocess"]["sparsify_dfs"]["after"]["m"] = g.count_edges();
			j["preprocess"]["sparsify_dfs"]["after"]["m_reduction"] = 1-1.0*g.count_edges()/cn;
		}
		Graph original_g = {1};
		std::unique_ptr<contract_tree_recovery_data> contracted1;
		std::unique_ptr<contract_tree_recovery_data> contracted2;
		if(contract_trees) {
			original_g = g;
			j["preprocess"]["contract_tree"]["before"]["n"] = g.n;
			j["preprocess"]["contract_tree"]["before"]["m"] = g.count_edges();
			int cn = g.n; 
			int cm = g.count_edges();
			stopwatch contract_time = {};
			contracted1 = contract_tree(g);
			swap(contracted1->g.edge_in, contracted1->g.edge_out);
			contracted2 = contract_tree(contracted1->g);
			swap(contracted2->g.edge_in, contracted2->g.edge_out);
			log_time(contract_time.total(), j["contract_tree"]["time"]);
			j["preprocess"]["contract_tree"]["after"]["n"] = contracted2->g.n;
			j["preprocess"]["contract_tree"]["after"]["m"] = contracted2->g.count_edges();
			g = contracted2->g;
			j["preprocess"]["contract_tree"]["after"]["n_reduction"] = 1-1.0*g.n/cn;
			j["preprocess"]["contract_tree"]["after"]["m_reduction"] = 1-1.0*g.count_edges()/cm;
		}
		////////
		auto tot = sw_s.lap();
		std::unique_ptr<Flowgraph<Edge::Minflow>> rg;
		if(solver_s == "pflowk2") {
			sw_s.lap();
			rg = pflowk2(g);
			auto tot = sw_s.lap();
			log_time(tot, j["solver"]["time"]);
		} else if(solver_s == "pflowk3") {
			sw_s.lap();
			rg = pflowk3(g);
			auto tot = sw_s.lap();
			log_time(tot, j["solver"]["time"]);
		} else {
			auto solver1 = std::find_if(mifsol.begin(), mifsol.end(), [&solver_s](auto u){return solver_s == u.second;});
			auto solver2 = std::find_if(mafsol.begin(), mafsol.end(), [&solver_s](auto u){return solver_s == u.second;});
			if(solver1 == mifsol.end() && solver2 == mafsol.end()) {
				std::cerr << "Unknown flow slover" << std::endl;
				exit(1);
			}
			tot = sw_s.lap();
			if(reduction_s == "naive") {
				rg = naive_minflow_reduction(g);
			} else if(reduction_s == "greedy") {
				rg = greedy_minflow_reduction(g);
			} else if(reduction_s == "greedy_sparsified") {
				rg = greedy_minflow_reduction_sparsified(g);
			} else {
				std::cerr << "Unknown minimum flow reduction" << std::endl;
				exit(1);
			}
			tot = sw_s.lap();
			log_time(tot, j["reduction"]["time"]);
			Flowgraph<Edge::Minflow> fgc(*rg);
			auto reduction_cover = minflow_reduction_path_recover_faster(fgc);
			j["reduction"]["cover"]["width"] = reduction_cover.size();
			sw_s.lap();
			if(solver1 != mifsol.end()) {
				solver1->first(*rg);
			} else {
				minflow_maxflow_reduction(*rg, solver2->first);
			}
			auto tot = sw_s.lap();
			log_time(tot, j["solver"]["time"]);
		}
		if(cover_decomposition == "pathcover") {
			tot = sw_s.lap();
			auto cover = minflow_reduction_path_recover_faster(*rg);
			tot = sw_s.lap();
			log_time(tot, j["path_recover"]["time"]);
			assert(is_valid_cover(cover, g));
			if(contract_trees) {
				sw_s.lap();
				swap(contracted2->g.edge_in, contracted2->g.edge_out);
				for(auto &p:cover)
					std::reverse(p.begin(), p.end());
				auto rp1 = recover_contract_pathcover(cover, *contracted2);
				for(auto &p:rp1)
					std::reverse(p.begin(), p.end());
				swap(contracted1->g.edge_in, contracted1->g.edge_out);
				cover = recover_contract_pathcover(rp1, *contracted1);
				log_time(sw_s.lap(), j["preprocess"]["contract_tree_recover"]["time"]);
				assert(is_valid_cover(cover, original_g));
			}
			j["cover"]["width"] = cover.size();
			int cover_len = 0;
			for(auto &u:cover)
				cover_len += u.size();
			j["cover"]["size"] = cover_len;
			j["memory"] = mem_peak();
			j["result"] = "ok";
			std::cout << j.dump() << std::endl;
			if(output_path.size() > 0) {
				std::ofstream out(output_path);
				if(!out.good()) {
					std::cerr << "Failed to open " << output_path << " for writing\n";
					exit(1);
				}
				out << cover.size() << "\n";
				for(auto &path:cover) {
					for(auto &u:path) {
						out << u << " ";
					}
					out << "\n";
				}
				out.close();
			}
		} else if(cover_decomposition == "chaincover_naive" || cover_decomposition == "chaincover_dict") {
			tot = sw_s.lap();
			chain_cover ccover;
			if(cover_decomposition == "chaincover_naive") {
				auto cover = minflow_reduction_path_recover_faster(*rg);
				ccover = *naive_chaincover_from_pathcover(cover, g);
			} else {
				ccover = *minflow_reduction_cc_naive(*rg);
			}
			tot = sw_s.lap();
			log_time(tot, j["cover_recover"]["time"]);
			assert(valid_cover(cover, g));
			assert(valid_chaincover(ccover, g));
			if(contract_trees) {
				std::cerr << "contract tree not supported with chaincover recovery";
				exit(1);
			}
			j["cover"]["width"] = ccover.size();
			int cover_len = 0;
			for(auto &u:ccover)
				cover_len += u.size();
			j["cover"]["size"] = cover_len;
			j["memory"] = mem_peak();
			j["result"] = "ok";
			std::cout << j.dump() << std::endl;
			if(output_path.size() > 0) {
				std::ofstream out(output_path);
				if(!out.good()) {
					std::cerr << "Failed to open " << output_path << " for writing\n";
					exit(1);
				}
				out << ccover.size() << "\n";
				for(auto &path:ccover) {
					for(auto &u:path) {
						out << u << " ";
					}
					out << "\n";
				}
				out.close();
			}
		} else {
			std::cerr << "Unknown cover decomposition: " << cover_decomposition << "\n";
			exit(1);
		}
		exit(0);
	} else {
		int status = 0;
		waitpid(pid, &status, 0); 
		auto tu = sw_c.total();
		if(WEXITSTATUS(status) == 0 && WIFEXITED(status)) { // normal
			return;
		}
		if(WEXITSTATUS(status) == 124) { // timeout
			j["result"] = "timeout";
		} else {
			j["result"] = "error";
			j["error"]["code"] = WEXITSTATUS(status);
		}
		log_time(tu, j["time"]);
		j["time"]["real"] = tu.real;
		std::cout << j.dump() << std::endl;
	}
}

Graph* read_graph2(std::string filename) {
	std::ifstream input(filename);
	if(!input.good()) {
		std::cerr << "Failed to open " << filename << " for reading\n";
	}
	int n, m;
	input >> n >> m;
	Graph *g = new Graph(n);
	for(int i=0; i<m; i++) {
		int a, b;
		input >> a >> b;
		g->add_edge(a, b);
	}
	input.close();
	return g;
}

std::vector<std::string> solvers() {
	std::vector<std::string> ret;
	ret.push_back("pflowk2");
	ret.push_back("pflowk3");
	for(auto &u:mifsol) {
		ret.push_back(u.second);
	}
	for(auto &u:mafsol) {
		ret.push_back(u.second);
	}
	return ret;
}

int main(int argc, char** args) {
	struct rlimit rlm;
	int ret = getrlimit(RLIMIT_STACK, &rlm);
	if(ret != 0)
		exit(1);
	if(rlm.rlim_cur != RLIM_INFINITY) {
		std::cerr << "Warning RLIMIT_STACK not set to infinity, but: " << rlm.rlim_cur << std::endl;
	}
	unsigned long mem_limit_bytes = std::numeric_limits<unsigned long>::max();
	int timeout_sec = std::numeric_limits<int>::max();
	bool contract_graph = false;
	std::string filepath = "";
	std::string reduction = "";
	std::string solver = "";
	std::string output_path = "";
	std::string cover_decomposition = "";
	int seed = 1337;
	bool sparsify_dfs_f = false;
	bool get_transitive_reduction = false;
	bool get_transitive_closure = false;
	int N, M, K;
	CLI::App app{"Example usage: ./exp -f random_dag -N 20 -M 25 -s pflowk2 --output_cover out_cover"};
	app.add_option("-f",filepath,"Either a path to a file, which is of format\n\t<# of nodes> <# of edges>\n\ta line for each edge a->b of the form <a b>\nOR\nrandom_dag (N and M must be provided)\nOR\nrandom_x_chain (N, M, K must be provided)")->required();
	app.add_option("-m",mem_limit_bytes,"memory limit in MB")->default_val(std::numeric_limits<unsigned long>::max());
	app.add_option("-t",timeout_sec,"time in sec")->default_val(std::numeric_limits<int>::max());
	app.add_option("-r",reduction,"reduction to use")
		->check(CLI::IsMember({"naive", "greedy", "greedy_sparsified"}));
	app.add_option("-s",solver,"solver to use")->required()
		->check(CLI::IsMember(solvers()));
	app.add_option("--contract_graph",contract_graph, "contract graph in out tree")->default_val(false);
	app.add_option("--sparsify_dfs",sparsify_dfs_f,"sparsify graph beforehand using the dfs sparsification")->default_val(false);
	app.add_option("--seed", seed, "seed to use when generating graphs")->default_val(1337);
	app.add_option("--transitive_reduction", get_transitive_reduction, "Use transitive reduction of the graph")->default_val(false);
	app.add_option("--transitive_closure", get_transitive_closure, "Use transitive closure of the graph")->default_val(false);
	app.add_option("-N",N,"N parameter for the graph generation");
	app.add_option("-M",M,"M parameter for the graph generation");
	app.add_option("-K",K,"K parameter for the graph generation");
	app.add_option("--cover_decomposition",cover_decomposition,"pathcover decomposes a path cover, chaincover_naive and dict decompose a chain cover by different methods")->check(CLI::IsMember({"pathcover", "chaincover_naive","chaincover_dict"}))->default_val("pathcover");
	app.add_option("--output_cover",output_path,"Optionally print the cover to the given path in format:\n\t<width>\t\none line for each path/chain with space separated nodes")->default_val("");
	CLI11_PARSE(app, argc, args);
	if(mem_limit_bytes != std::numeric_limits<unsigned long>::max()) {
		mem_limit_bytes *= 1024*1024;
	}
	
	std::unique_ptr<Graph> g;
	nlohmann::json j;
	if(filepath == "random_dag") {
		g = random_dag(N, M, seed);
	} else if(filepath == "funnel_gen") {
		g = funnel_gen(N);
	} else if(filepath == "random_x_partite") {
		g = random_x_partite(log(N), N, M, seed); // lower bound ?
	} else if(filepath == "random_x_chain") {
		g = random_x_chain(K, N, M, seed); // upper bound
		j["graph"]["K_PARAM"] = K;
	} else {
		g = std::make_unique<Graph>(*read_graph2(filepath));
	}
	if(get_transitive_closure) {
		g = transitive_closure(*g);
	}
	if(get_transitive_reduction) {
		g = transitive_reduction(*g);
	}
	j["graph"]["name"] = filepath;
	j["graph"]["N_PARAM"] = N;
	j["graph"]["M_PARAM"] = M;
	run_one(*g, timeout_sec, mem_limit_bytes, reduction, solver, sparsify_dfs_f, contract_graph, output_path, cover_decomposition, j);
}
