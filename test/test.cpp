#include <mpc/naive.h>
#include <mpc/graph.h>
#include <mpc/lemon.h>
#include <mpc/cc.cpp>
#include <mpc/reach.h>
#include <mpc/antichain.cpp>
#include <mpc/pflow.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include <mpc/preprocessing.h>

struct test_graph {
	Graph *g;
	int width;
};

Graph* read_graph2(std::string filename) {
	std::ifstream input(filename);
	assert(input.good());
	int width;
	input >> width;
	int n, m;
	input >> n >> m;
	Graph *g = new Graph(n);
	for(int i=0; i<m; i++) {
		int a, b;
		input >> a >> b;
		g->add_edge(a, b);
	}
	return g;
}

test_graph read_graph(std::string filename) {
	std::ifstream input(filename);
	assert(input.good());
	int width;
	input >> width;
	int n, m;
	input >> n >> m;
	Graph *g = new Graph(n);
	for(int i=0; i<m; i++) {
		int a, b;
		input >> a >> b;
		g->add_edge(a, b);
	}
	return {g, width};
}

void reach_test(Graph &g, path_cover &pc) {
	auto ri = reachability_idx(g, pc);
	for(int i=1; i<g.n; i++) {
		for(int j=1; j<=g.n; j++) {
			ASSERT_TRUE(naive_reach(i, j, g) == ri.reaches(i, j)) << i << " " << j;
		}
	}
}

void test_dfs_sparsify(Graph &g) {
	auto dg = sparsify_dfs_elegant(g);
	auto g1 = pflowk2(g);
	auto g2 = pflowk2(*dg);
	auto pc1 = minflow_reduction_path_recover_faster(*g1);
	auto pc2 = minflow_reduction_path_recover_faster(*g2);
	ASSERT_TRUE(pc1.size() == pc2.size() && is_valid_cover(pc2, g));
	auto r1 = reachability_idx(g, pc1);
	auto r2 = reachability_idx(*dg, pc1);
	for(int i=1; i<=g.n; i++)
		for(int j=1; j<=g.n; j++)
			ASSERT_TRUE(r1.reaches(i,j) == r2.reaches(i,j));
}

void test_all(Graph &g) {
	std::vector<std::function<void(Flowgraph<Edge::Minflow>&)>> mifsol = {lemon_cs, lemon_cs_minlen, lemon_ns, lemon_ns_minlen, lemon_cc, lemon_cc_minlen, lemon_caps, lemon_caps_minlen, naive_minflow_solve};
	std::vector<std::function<void(Flowgraph<Edge::Maxflow>&)>> mafsol = {lemon_preflow, maxflow_solve_edmonds_karp, maxflow_solve_edmonds_karp_DMOD};

	std::vector<path_cover> pc;
	int prev_width = -1;
	auto pflowc = pflowk2(g);
	auto pflowpc = minflow_reduction_path_recover_faster(*pflowc);
	EXPECT_TRUE(is_valid_cover(pflowpc, g));
	pc.push_back(pflowpc);
	auto test = [&mafsol, &mifsol, &pc, &g](Flowgraph<Edge::Minflow> &rgo) {
		for(auto sol:mifsol) {
			auto rg = std::make_unique<Flowgraph<Edge::Minflow>>(rgo);
			sol(*rg);
			//
			auto mac = maxantichain_from_minflow(*rg);
			//
			auto cover = minflow_reduction_path_recover_faster(*rg);
			EXPECT_TRUE(is_valid_cover(cover, g));
			pc.push_back(cover);

			EXPECT_TRUE(mac.size() == cover.size()) << "MaxAntichain size" << mac.size() <<  " mpc width: " << cover.size();
			EXPECT_TRUE(is_antichain(mac, g));

			auto chain_cover1 = naive_chaincover_from_pathcover(cover, g);
			auto rg2 = std::make_unique<Flowgraph<Edge::Minflow>>(rgo);
			sol(*rg2);
			auto chain_cover2 = minflow_reduction_cc_fast(*rg2);
			EXPECT_TRUE(valid_chaincover(*chain_cover1, g));
			EXPECT_TRUE(valid_chaincover(*chain_cover2, g));
			EXPECT_TRUE(chain_cover1->size() == chain_cover2->size());
			reach_test(g, cover); 
		}
		for(auto sol:mafsol) {
			auto rg = std::make_unique<Flowgraph<Edge::Minflow>>(rgo);
			minflow_maxflow_reduction(*rg, sol);
			auto cover = minflow_reduction_path_recover_faster(*rg);
			EXPECT_TRUE(is_valid_cover(cover, g));
			pc.push_back(cover);
		}
	};
	auto r1 = naive_minflow_reduction(g);
	test(*r1);
	auto r2 = greedy_minflow_reduction(g);
	test(*r2);
	auto r4 = greedy_minflow_reduction_sparsified(g);
	test(*r4);
	for(int i=1; i<pc.size(); i++) {
		ASSERT_TRUE(pc[i-1].size() == pc[i].size()) << pc[i-1].size() << " " << pc[i].size();
	}
}

void random_dag_test_all(int seed) {
	auto g = random_dag(18, 15, 109624);
	test_all(*g);


	for(int m=0; m<=1000; m+=10) {
		auto g = random_dag(100, m, seed);
		//g->print();
		test_all(*g);
	}
	for(int m=1000; m<=4000; m+=1000) {
		auto g = random_dag(100, m, seed);
		test_all(*g);
	}
}

void random_partite_test_all(int seed) {
	for(int m=0; m<=1000; m+=10) {
		for(int x=2; x<10; x++) {
			auto g = random_x_partite(x, 10, m, seed);
			test_all(*g);
		}
	}
	for(int m=1000; m<=4000; m+=1000) {
		for(int x=2; x<10; x++) {
			auto g = random_x_partite(x, 10, m, seed);
			test_all(*g);
		}
	}
}

void random_chain_test_all(int seed) {
	for(int k=2; k<10; k++) {
		int N = 10;
		for(int m=0; m<5; m++) {
			auto g = random_x_chain(k, N, N*m, seed);
			test_all(*g);
		}
	}
}

void random_binarytree_test_all(int seed) {
	for(int i=1; i<=8; i++) {
		auto g = binary_tree(i, false);
		test_all(*g);
		g = binary_tree(i, true);
		test_all(*g);
	}
}

class tc1 :
    public testing::TestWithParam<int> {
};
TEST_P(tc1, random_dag_test_all) {
	random_dag_test_all(GetParam());
}
TEST_P(tc1, random_binarytree_test_all) {
	random_binarytree_test_all(GetParam());
}
TEST_P(tc1, random_partite_test_all) {
	random_partite_test_all(GetParam());
}
TEST_P(tc1, random_chain_test_all) {
	random_chain_test_all(GetParam());
}

TEST_P(tc1, dfs_sparsify) {
	for(int i=1; i<=8; i++) {
		auto g = binary_tree(i, false);
		test_dfs_sparsify(*g);
		g = binary_tree(i, true);
		test_dfs_sparsify(*g);
	}
	for(int m=0; m<=(100*99/2); m+=(m < 500 ? 40 : 400)) {
		auto g = random_dag(100, m, GetParam());
		test_all(*g);
		test_dfs_sparsify(*g);
	}
}
INSTANTIATE_TEST_SUITE_P(tc, tc1, ::testing::Range(1, 10));

class tc2 :
    public testing::TestWithParam<int> {
};

TEST_P(tc2, naive_minflow_reduction_width_test) {
	auto g = read_graph("data/test" + std::to_string(GetParam()));
	auto red_fg = naive_minflow_reduction(*g.g);
	naive_minflow_solve(*red_fg);
	auto cover = minflow_reduction_path_recover(*red_fg);
	ASSERT_TRUE(cover.size() == g.width);
	test_all(*g.g);
}

INSTANTIATE_TEST_SUITE_P(tc, tc2, ::testing::Range(1, 5));
